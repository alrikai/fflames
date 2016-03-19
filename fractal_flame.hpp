/* fractal_flame.hpp -- part of the fractal flames implementation 
 *
 * Copyright (C) 2015 Alrik Firl 
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef FF_FRACTAL_FLAME_HPP
#define FF_FRACTAL_FLAME_HPP

#include "ifs.hpp"
#include "ifs_types.hpp"
#include "util/ff_utils.hpp"

#include <opencv2/opencv.hpp>

#include <random>
#include <iostream>
#include <utility>
#include <map>
#include <chrono>
#include <thread>
#include <mutex>

template <typename data_t, typename pixel_t>
class fflame_data
{
public:    
    //using hist_t = std::array<std::array<histogram_info<pixel_t>, fflame_constants::imwidth>, fflame_constants::imheight>;
    using hist_t = std::vector<histogram_info<pixel_t>>;

    fflame_data(int h_height = fflame_constants::imheight, int h_width = fflame_constants::imwidth)
        : hist_height(h_height), hist_width(h_width), fflame_hist(hist_t(hist_height * hist_width, histogram_info<pixel_t>())) 
    {} 

    //NOTE: no copy-constructor, since std::mutex is move-only
    fflame_data(const fflame_data<data_t, pixel_t>&) = delete;
    fflame_data<data_t, pixel_t>& operator= (const fflame_data<data_t, pixel_t>&) = delete;
  
    template <class hist_container_t>
    void apply_fflame_run(hist_container_t&& hist_data)
    {
        std::lock_guard<std::mutex> lk (histdata_mtx);

        for (int row = 0; row < hist_height; row++) {
            for (int col = 0; col < hist_width; col++) {
                int idx = row * hist_width + col;
                fflame_hist[idx].update(hist_data[idx]);    
            }
        }
    }
    
    void get_and_reset(hist_t& hist_data)
    {
        std::lock_guard<std::mutex> lk (histdata_mtx);

        std::copy(fflame_hist.begin(), fflame_hist.end(), hist_data.begin());
        std::for_each(fflame_hist.begin(), fflame_hist.end(), 
                [](histogram_info<pixel_t>& hdata)
                {
                    hdata.reset();
                });
    }

private:
    int hist_height;
    int hist_width;

    std::mutex histdata_mtx;
    hist_t fflame_hist;

    //cv::Mat_<histogram_inf>>pixel_t>> fflame_hist;
};

template <typename data_t, typename pixel_t>
void run_fflame(const affine_fcns::invoker<data_t>* const flamer, const int num_points, fflame_data<data_t, pixel_t>* fdata, fflame_util::fast_rand& f_rand)
{
    //NOTE: using an std::array is technically possible here, but if the image size gets large then we'll have stack issues
    //std::array<histogram_info<pixel_t>, fflame_constants::imheight * fflame_constants::imwidth> hist_data;
    std::vector<histogram_info<pixel_t>> hist_data (fflame_constants::imheight * fflame_constants::imwidth);

    //for generating the random point to use
    std::uniform_real_distribution<> dis(fflame_constants::min_pt, fflame_constants::max_pt);
    auto fflame_rngengine = fflame_util::get_engine();

    const data_t cosidx_rotate_factor = std::cos(2 * fflame_constants::PI);
    const data_t sinidx_rotate_factor = std::sin(2 * fflame_constants::PI);
    const auto fflame_point_range = fflame_constants::max_pt - fflame_constants::min_pt;

    for (int sample = 0; sample < num_points; ++sample)
    {
        flame_point<data_t> flame_pt (dis(fflame_rngengine), dis(fflame_rngengine));
        for (int i = 0; i < fflame_constants::max_iter; ++i)
        {
            //NOTE: if we have per-function probability weights, would have to use them here to select the function
            int fcn_idx = f_rand.xorshift128plus(flamer->fcn.size());
            flamer->invoke(fcn_idx, flame_pt);
         
            //don't store the first 20 iterations, as they'll apparently be too far from the solution
            if(i > 20)
            {
                //store the intermediate result to the histogram datastructures
                data_t rotated_x = flame_pt.x * cosidx_rotate_factor - flame_pt.y * sinidx_rotate_factor;
                data_t rotated_y = flame_pt.x * sinidx_rotate_factor + flame_pt.y * cosidx_rotate_factor;
                long col_idx = std::lround(fflame_constants::imwidth - ((fflame_constants::max_pt - rotated_x) / fflame_point_range) * fflame_constants::imwidth);
                long row_idx = std::lround(fflame_constants::imheight - ((fflame_constants::max_pt - rotated_y) / fflame_point_range) * fflame_constants::imheight);

                //if the point (at image dimensions) is within bounds, add it to the histogram 
                if((col_idx >= 0 && col_idx < fflame_constants::imwidth) && (row_idx >= 0 && row_idx < fflame_constants::imheight))
                {
                    pixel_t color (flame_pt.color[0], flame_pt.color[1], flame_pt.color[2]);
                    hist_data[row_idx*fflame_constants::imwidth+col_idx].update(color, 1);
                }
            }
        }
    }

    //apply the data gathered from this run to the overall histograms
    fdata->apply_fflame_run(std::move(hist_data));
}



#if 0
template <typename data_t, typename pixel_t = cv::Vec<data_t, 3>>
void generate_fractal_flame(affine_fcns::invoker<data_t>& flamer, cv::Mat_<pixel_t>& image, const int num_workers = std::thread::hardware_concurrency())
{ 
    //start timing (just for informative purposes)
    auto start_time = std::chrono::high_resolution_clock::now();
 
    //for seeding the flame thread rng 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> fast_dis(0, std::numeric_limits<uint64_t>::max());

    std::vector<std::unique_ptr<flame_thread>> fflame_tasks;
    //spawn the threads...
    for (int th_count = 0; th_count < num_workers; ++th_count) { 
        fflame_tasks.emplace_back(new flame_thread(fast_dis(gen), fast_dis(gen)));
        bool flame_run = fflame_tasks.at(th_count)->do_flame(&run_fflame<data_t, pixel_t>, flamer, fflame_constants::num_samples/num_workers);
        //in light of the fact that this really shouldn't happen, should I be throwing an exception here?
        if(!flame_run) {
            std::cout << "NOTE: flame task failed" << std::endl;
        }
    }
        
    //...join the threads
    for (int th_count = 0; th_count < num_workers; ++th_count) {
        fflame_tasks.at(th_count)->finish_flame();
    }

    /*
     *  @this point -- we will have finished all of the fractal flame generation steps and 
     *  created the histograms above. From here on down it'll be single threaded (for now),
     *  so it might look weird, but we'll just access the histograms and image in fflame_data 
     *  without any manner of synchronization, since we know all the above work is done and
     *  it's all single-threaded execution from here on down
     */
    std::cout << "Done with generation, moving to image generation..." << std::endl;
    image = cv::Mat_<pixel_t>::zeros(fflame_constants::imheight, fflame_constants::imwidth);
    render_fractal_flame<data_t, pixel_t>(image);

    auto current_time = std::chrono::high_resolution_clock::now();
    double time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
    std::cout << "Took " << time_elapsed << " ms" << std::endl;
}
#endif




#endif
