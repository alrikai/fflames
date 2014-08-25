
#include "ifs.hpp"
#include "ifs_types.hpp"
//#include "thread_pool.hpp"

#include <opencv2/opencv.hpp>

#include <random>
#include <iostream>
#include <utility>
#include <chrono>

#include <future> 


template <typename data_t, typename pixel_t>
class fflame_data
{
public:    
    fflame_data()
    {
        freq_hist = cv::Mat_<data_t>::zeros(fflame_constants::imheight, fflame_constants::imwidth);
        color_hist = cv::Mat_<pixel_t>::zeros(fflame_constants::imheight, fflame_constants::imwidth);
    }

    //NOTE: no copy-constructor, since std::mutex is move-only
    fflame_data(const fflame_data<data_t, pixel_t>&) = delete;
    fflame_data<data_t, pixel_t>& operator= (const fflame_data<data_t, pixel_t>&) = delete;

    void apply_fflame_run(std::map<std::pair<size_t, size_t>, histogram_info<pixel_t>>&& hist_data)
    {
        std::lock_guard<std::mutex> lk (histdata_mtx);
        auto hist_it = hist_data.begin();
        while(hist_it != hist_data.end())
        {
            const int row_idx = hist_it->first.first;
            const int col_idx = hist_it->first.second;
    
            freq_hist(row_idx, col_idx) += hist_it->second.frequency_count;
            color_hist(row_idx, col_idx) = (color_hist(row_idx, col_idx) + hist_it->second.color) / 2.0;

            hist_it++;
        }
    }

//TODO: this part is likely temporary -- since we are using single-threaded image generation, as it's a negligible amount of time to run
//private:
    mutable std::mutex histdata_mtx;
	cv::Mat_<data_t> freq_hist;
	cv::Mat_<pixel_t> color_hist;
};


/////////////////////////////////////////////////////////////////////////////
//NOTE: need to re-arrange everything here. Am making this a global object since it has 
//to be used by the threads in run_fflame, but fflame_data is non-copyable on account of
//the mutex data member, hence I can't pass references to the object around. Should I just
//pass a pointer around instead?
/////////////////////////////////////////////////////////////////////////////
template <typename data_t, typename pixel_t>
struct fflame_histdata
{
    static fflame_data<data_t, pixel_t> fdata;
};
template <typename data_t, typename pixel_t>
fflame_data<data_t, pixel_t> fflame_histdata<data_t, pixel_t>::fdata;
/////////////////////////////////////////////////////////////////////////////



template <typename data_t, typename pixel_t>
void run_fflame(const affine_fcns::invoker<data_t>& flamer, const int num_points)
{
    std::map<std::pair<size_t, size_t>, histogram_info<pixel_t>> hist_data;
    //for selecting the flame function to use
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> fast_dis(0, std::numeric_limits<uint64_t>::max());
    fflame_randutil::fast_rand f_rand(fast_dis(gen), fast_dis(gen));

    //for generating the random point to use
    std::uniform_real_distribution<> dis(fflame_constants::min_pt, fflame_constants::max_pt);

    for (int sample = 0; sample < num_points; ++sample)
    {
        flame_point<data_t> flame_pt (dis(fflame_randutil::get_engine()), dis(fflame_randutil::get_engine()));
        for (int i = 0; i < fflame_constants::max_iter; ++i)
        {
            //NOTE: if we have per-function probability weights, would have to use them here to select the function
            int fcn_idx = f_rand.xorshift128plus(affine_fcns::invoker<data_t>::fcn.size()); //(math_constants::get_engine()); 
            flamer.invoke(fcn_idx, flame_pt);
         
            //don't store the first 20 iterations, as they'll apparently be too far from the solution
            if(i > 20)
            {
                //store the intermediate result to the histogram datastructures
                
                double rotated_x = flame_pt.x * std::cos(2 * fflame_constants::PI) - flame_pt.y * std::sin(2 * fflame_constants::PI);
                double rotated_y = flame_pt.x * std::sin(2 * fflame_constants::PI) + flame_pt.y * std::cos(2 * fflame_constants::PI);
                
                long col_idx = std::lround(fflame_constants::imwidth - ((fflame_constants::max_pt - rotated_x) / (fflame_constants::max_pt - fflame_constants::min_pt)) * fflame_constants::imwidth);
                long row_idx = std::lround(fflame_constants::imheight - ((fflame_constants::max_pt - rotated_y) / (fflame_constants::max_pt - fflame_constants::min_pt)) * fflame_constants::imheight);

                //if the point (at image dimensions) is within bounds, add 
                if((col_idx >= 0 && col_idx < fflame_constants::imwidth) && (row_idx >= 0 && row_idx < fflame_constants::imheight))
                {
                    auto hist_idx = std::make_pair(row_idx, col_idx);
                    pixel_t color (flame_pt.color[0], flame_pt.color[1], flame_pt.color[2]);
                    
                    //update the map -- add to an existing entry or generate a new one
                    auto bin_entry = hist_data.find(hist_idx); 
                    if(bin_entry != hist_data.end())
                        bin_entry->second.update(color, 1);
                    else
                        hist_data.insert(std::make_pair(hist_idx,histogram_info<pixel_t>(color, 1)));
                }
            }
        }
    }

    //apply the data gathered from this run to the overall histograms
    fflame_histdata<data_t, pixel_t>::fdata.apply_fflame_run(std::move(hist_data));
}

int main()
{
    using data_t = double;
    using pixel_t = cv::Vec<data_t, 3>;

    affine_fcns::invoker<data_t> flamer;
    flamer.randomize_parameters(-2, 2); 

    const int num_workers = std::thread::hardware_concurrency();
  
    //start timing (just for informative purposes)
    auto start_time = std::chrono::high_resolution_clock::now();
  
    std::vector<std::thread> fflame_tasks;
    //spawn the threads...
    for (int th_count = 0; th_count < num_workers; ++th_count)
        fflame_tasks.push_back(std::thread(&run_fflame<data_t, pixel_t>, flamer, fflame_constants::num_samples/num_workers));    
        
    //...join the threads
    for (int th_count = 0; th_count < num_workers; ++th_count)
        fflame_tasks[th_count].join();

/*
    for (int sample = 0; sample < fflame_constants::num_samples; sample+=num_workers)
    {
        flame_point<data_t> flame_pt (dis(math_constants::get_engine()), dis(math_constants::get_engine()));
        run_fflame<data_t, pixel_t>(flamer, std::move(flame_pt));
    }
*/

    /*
     *  @this point -- we will have finished all of the fractal flame generation steps and 
     *  created the histograms above. From here on down it'll be single threaded (for now),
     *  so it might look weird, but we'll just access the histograms and image in fflame_data 
     *  without any manner of synchronization, since we are assuming all the above work is 
     *  done and single-threaded execution from here on down
     */
    std::cout << "Done with generation, moving to image generation..." << std::endl;
   
    cv::Mat_<pixel_t> image = cv::Mat_<data_t>::zeros(fflame_constants::imheight, fflame_constants::imwidth);
    auto& freq_hist = fflame_histdata<data_t, pixel_t>::fdata.freq_hist;
    auto& color_hist = fflame_histdata<data_t, pixel_t>::fdata.color_hist;

    const data_t freq_max_log = std::log10(*std::max_element(freq_hist.begin(), freq_hist.end()));

    //NOTE: we assume the histogram dimensions are multiples of the image dimensions
    const int rowpx_factor = fflame_constants::hist_height/fflame_constants::imheight;
    const int colpx_factor = fflame_constants::hist_width/fflame_constants::imwidth;
    //make the final output image
    for (int im_row = 0; im_row < fflame_constants::imheight; ++im_row)
    {
        int hist_row = im_row*rowpx_factor;     
        for (int im_col = 0; im_col < fflame_constants::imwidth; ++im_col)
        {
            data_t freq_avg = 0;
            pixel_t color_avg = 0;
            int hist_col = im_col*colpx_factor;    
            for (int supersample_row = 0; supersample_row < rowpx_factor; ++supersample_row)
            {
                for (int supersample_col = 0; supersample_col < colpx_factor; ++supersample_col)
                {
                    freq_avg += freq_hist(hist_row + supersample_row, hist_col + supersample_col);
                    color_avg += color_hist(hist_row + supersample_row, hist_col + supersample_col);
                }
            }

            freq_avg /= rowpx_factor*colpx_factor;
            color_avg /= rowpx_factor*colpx_factor;

            data_t alpha = std::log10(freq_avg)/freq_max_log;
            image(im_row, im_col) = 255 * color_avg * std::pow(alpha, fflame_constants::gamma_factor); 
        }
    }


    auto current_time = std::chrono::high_resolution_clock::now();
    double time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
    std::cout << "Took " << time_elapsed << " ms" << std::endl;


    cv::imwrite("image.png", image);
}


