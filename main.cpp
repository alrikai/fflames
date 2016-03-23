/* main.cpp -- part of the fractal flames implementation 
 *
 * Copyright (C) 2015 Alrik Firl 
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */



#include "fflame_generator.hpp"
#include "fractal_flame.hpp"

#include <opencv2/opencv.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <memory>
#include <stdexcept>

namespace bfs = boost::filesystem;

//--------------------------------------------------------------------------------------------------------
//for debugging with gdb
double get_mat_pixel(cv::Mat_<cv::Vec<double,3>>& image, int row, int col, int channel)
{
    auto px = image(row, col);
    return px(channel);
}
//--------------------------------------------------------------------------------------------------------

//generates num_frames frames in interpolating between the lhs and rhs frames
template <typename pixel_t>
int interpolate_frames(const cv::Mat_<pixel_t>& lhs_img, const cv::Mat_<pixel_t>& rhs_img, const int num_frames, const bfs::path out_filebasepath, int frame_idx)
{
    auto output_fpath = out_filebasepath;
    output_fpath /= std::to_string(frame_idx++) + ".png";
    const std::string out_filepath = output_fpath.native();

    {
    //NOTE: we'll 'lose' the last frame generated in this manner?
    auto output_keyframe_fpath = out_filebasepath;
    output_keyframe_fpath /= "baseframe_" + std::to_string(frame_idx++) + ".png";
    const std::string out_keyframe_filepath = output_keyframe_fpath.native();
    cv::imwrite(out_keyframe_filepath, lhs_img);
    }

    //std::cout << "saving " << out_filepath << std::endl;
    cv::imwrite(out_filepath, lhs_img);

    cv::Mat_<pixel_t> pixel_diffimg = (rhs_img - lhs_img) / static_cast<double>(num_frames);    
    //note: this is just a reference to lhs_img, not a clone
    cv::Mat_<pixel_t> working_frame = lhs_img;
    for (int i = 0; i < num_frames; ++i) {
        working_frame += pixel_diffimg;
        output_fpath = out_filebasepath;
        output_fpath /= std::to_string(frame_idx++) + ".png";
        const std::string out_filepath = output_fpath.native();
        cv::imwrite(out_filepath, working_frame);
        //std::cout << "saving " << out_filepath << std::endl;
    }

    //at this point, working_frame should be 1 step off from rhs_img
    return frame_idx;
}

template <typename pixel_t>
int interpolate_frames(const flame_frame<pixel_t>& lhs_frame, const flame_frame<pixel_t>& rhs_frame, const int num_frames, const bfs::path out_filebasepath, int frame_idx)
{
    cv::Mat_<pixel_t> lhs_img (lhs_frame.rows, lhs_frame.cols, lhs_frame.data);
    cv::Mat_<pixel_t> rhs_img (rhs_frame.rows, rhs_frame.cols, rhs_frame.data);
    return interpolate_frames (lhs_img, rhs_img, num_frames, out_filebasepath, frame_idx);
}



template <template <class> class frame_t, typename pixel_t, typename data_t>
void generate_fractal_flames (const bfs::path output_dir, const int num_working_variants, const int num_images, const bool do_interpolation = true)
{
    using flame_gen_t = fflame_generator<frame_t, data_t, pixel_t>;
    const int num_generator_threads = 4;
    auto bg_generator = std::unique_ptr<flame_gen_t> (new flame_gen_t(fflame_constants::imheight, fflame_constants::imwidth, num_working_variants, num_generator_threads));
    //pass in a queue to hold the finished flame frames
    //auto bg_framequeue = std::unique_ptr<EventQueue<frame_t<pixel_t>>>(new EventQueue<frame_t<pixel_t>>(num_images));

    typename flame_gen_t::template ts_queue_t <frame_t<pixel_t>> bg_framequeue;
    auto bg_framequeue_p = &bg_framequeue;
    bg_generator->register_framequeue(bg_framequeue_p);    
    bg_generator->start_generation();

    int output_index = 0;
    std::unique_ptr<frame_t<pixel_t>> prev_image;
    //frame_t<pixel_t> prev_image;
    for (int frame_idx = 0; frame_idx < num_images; ++frame_idx) {

        bool got_frame = false;
        std::unique_ptr<frame_t<pixel_t>> fflame_frame;
        while (!got_frame) {
            fflame_frame = bg_framequeue.pop(got_frame);
        }
            
        if(got_frame && fflame_frame) {
            //save the frames out, in whatever configuration is specified
            if(do_interpolation) {
                if(frame_idx > 0) {
                    output_index = interpolate_frames(*(prev_image.get()), *(fflame_frame.get()), 30, output_dir, output_index);
                }
                prev_image.swap(fflame_frame);
            } else {
                cv::Mat_<pixel_t> keyframe_img (fflame_frame->rows, fflame_frame->cols, fflame_frame->data);
                //NOTE: we'll 'lose' the last frame generated in this manner?
                auto output_keyframe_fpath = output_dir;
                output_keyframe_fpath /= "baseframe_" + std::to_string(frame_idx) + ".png";
                const std::string out_keyframe_filepath = output_keyframe_fpath.native();
                cv::imwrite(out_keyframe_filepath, keyframe_img);
            }
        }
    }

    bg_generator->stop_generation();
}

//------------------------------------------------------------------------------------------------------
template <typename pixel_t>
using frame_t = flame_frame<pixel_t>; //cv::Mat_<pixel_t>;

int main(int argc, char* argv[])
{
    std::vector<std::string> manual_variant_list; 

    namespace bpo = boost::program_options; 
    bpo::options_description bpo_desc("fractal flames options"); 
    bpo_desc.add_options() 
        ("help,h", "Print help message")
        ("output-path,o", bpo::value<std::string>(), "output flame-image path")
        ("num-variants,n", bpo::value<int>(), "number of working variants")
        ("total-frames,t", bpo::value<int>(), "total number of keyframes to generate")
        ("variant-list,v", bpo::value<std::vector<std::string>>(&manual_variant_list), "list of variants to use")
        ("do-interpolation,i", bpo::value<int>(), "flag for writing out interpolated frames between keyframes (1: yes, 0: no)");

    bpo::variables_map vm;
    bpo::store(bpo::command_line_parser(argc, argv).options(bpo_desc).run(), vm);
    bpo::notify(vm);

    std::string output_path;
    int num_working_variants, num_images, do_interpolation;

    if(vm.count("help")){
        std::cout << bpo_desc << std::endl;
        return 0;
    }

    if(vm.count("output-path")){
        output_path = vm["output-path"].as<std::string>();
    } else {
        std::cout << bpo_desc << std::endl;
        return 0;
    }

    if(vm.count("num-variants")){
        num_working_variants = vm["num-variants"].as<int>();
    } else {
        std::cout << bpo_desc << std::endl;
        return 0;
    }

    if(vm.count("total-frames")){
        num_images = vm["total-frames"].as<int>();
    } else {
        std::cout << bpo_desc << std::endl;
        return 0;
    }

    if(vm.count("do-interpolation")){
        do_interpolation = vm["do-interpolation"].as<int>();
    } else {
        std::cout << bpo_desc << std::endl;
        return 0;
    }

    //invoke via something like
    //./fflame_gen -o ff_frames/ff_v13 -n 10 -t 100 -i 0 -v fisheye -v eyefish -v swirl -v this -v is -v a -v drill
    //
    //... where the -v are the manually provided variants, given in sequential order
    //next we just have to figure out how to pass these in...
    std::cout << "Given Variant (" << manual_variant_list.size() << ") List: " << std::endl;
    for (auto var_id : manual_variant_list) {
        std::cout << var_id << std::endl;
    }

    //check the output directory path (if it doesn't exist, create it)
    bfs::path output_dir(output_path);
    if (!(bfs::exists(output_dir) && bfs::is_directory(output_dir))) {
        if (!boost::filesystem::create_directory(output_path)) {
            std::cout << "Invalid output file directory -- " << output_path << std::endl;
            throw std::invalid_argument("Invalid output file directory " + output_path);
        }
    }

    using data_t = double;
    using pixel_t = cv::Vec<data_t, 3>;  

    bool do_interpolation_flag = do_interpolation != 0;
    generate_fractal_flames<frame_t, pixel_t, data_t> (output_dir, num_working_variants, num_images, do_interpolation_flag);

    return 0;
}

