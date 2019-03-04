/* main.cpp -- part of the fractal flames implementation
 *
 * Copyright (C) 2015 Alrik Firl
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include "fflame_generator.hpp"
#include "fractal_flame.hpp"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <opencv2/opencv.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace bfs = boost::filesystem;

//--------------------------------------------------------------------------------------------------------
// for debugging with gdb
double get_mat_pixel(cv::Mat_<cv::Vec<double, 3>> &image, int row, int col,
                     int channel) {
  auto px = image(row, col);
  return px(channel);
}
//--------------------------------------------------------------------------------------------------------

// generates num_frames frames in interpolating between the lhs and rhs frames
template <typename pixel_t>
int interpolate_frames(const cv::Mat_<pixel_t> &lhs_img,
                       const cv::Mat_<pixel_t> &rhs_img, const int num_frames,
                       const bfs::path out_filebasepath, int frame_idx) {
  auto output_fpath = out_filebasepath;
  output_fpath /= std::to_string(frame_idx++) + ".png";
  const std::string out_filepath = output_fpath.native();

  {
    // NOTE: we'll 'lose' the last frame generated in this manner?
    auto output_keyframe_fpath = out_filebasepath;
    output_keyframe_fpath /=
        "baseframe_" + std::to_string(frame_idx++) + ".png";
    const std::string out_keyframe_filepath = output_keyframe_fpath.native();
    cv::imwrite(out_keyframe_filepath, lhs_img);
  }

  // std::cout << "saving " << out_filepath << std::endl;
  cv::imwrite(out_filepath, lhs_img);

  cv::Mat_<pixel_t> pixel_diffimg =
      (rhs_img - lhs_img) / static_cast<double>(num_frames);
  // note: this is just a reference to lhs_img, not a clone
  cv::Mat_<pixel_t> working_frame = lhs_img;
  for (int i = 0; i < num_frames; ++i) {
    working_frame += pixel_diffimg;
    output_fpath = out_filebasepath;
    output_fpath /= std::to_string(frame_idx++) + ".png";
    const std::string out_filepath = output_fpath.native();
    cv::imwrite(out_filepath, working_frame);
    // std::cout << "saving " << out_filepath << std::endl;
  }

  // at this point, working_frame should be 1 step off from rhs_img
  return frame_idx;
}

template <typename pixel_t>
int interpolate_frames(const flame_frame<pixel_t> &lhs_frame,
                       const flame_frame<pixel_t> &rhs_frame,
                       const int num_frames, const bfs::path out_filebasepath,
                       int frame_idx) {
  cv::Mat_<pixel_t> lhs_img(lhs_frame.rows, lhs_frame.cols, lhs_frame.data);
  cv::Mat_<pixel_t> rhs_img(rhs_frame.rows, rhs_frame.cols, rhs_frame.data);
  return interpolate_frames(lhs_img, rhs_img, num_frames, out_filebasepath,
                            frame_idx);
}

template <class framequeue_t, template <class> class frame_t, typename pixel_t>
void run_frame_generation(framequeue_t &bg_framequeue,
                          const bfs::path output_dir, const int num_images,
                          bool do_interpolation) {
  int output_index = 0;
  std::unique_ptr<frame_t<pixel_t>> prev_image;
  // frame_t<pixel_t> prev_image;
  for (int frame_idx = 0; frame_idx < num_images; ++frame_idx) {

    bool got_frame = false;
    std::unique_ptr<frame_t<pixel_t>> fflame_frame;
    while (!got_frame) {
      fflame_frame = bg_framequeue.pop(got_frame);
    }

    if (got_frame && fflame_frame) {
      // save the frames out, in whatever configuration is specified
      if (do_interpolation) {
        if (frame_idx > 0) {
          output_index =
              interpolate_frames(*(prev_image.get()), *(fflame_frame.get()), 30,
                                 output_dir, output_index);
        }
        prev_image.swap(fflame_frame);
      } else {
        /*
        //Q: how do we normalize? We could get the maximum R, G, and B hcannels?
        pixel_t min_pxval;
        pixel_t max_pxval;
        std::for_each(fflame_frame->data, fflame_frame->data +
        fflame_frame->rows * fflame_frame->cols,
                    [&min_pxval, &max_pxval] (const pixel_t& pxval) {
                        for (int k = 0; k < pixel_t::channels; k++) {
                            if(pxval[k] > max_pxval[k]) {
                                max_pxval[k] = pxval[k];
                            }
                            if(pxval[k] < min_pxval[k]) {
                                min_pxval[k] = pxval[k];
                            }
                        }
                    });

        std::cout << "Per-channel MIN: {" << min_pxval[0] << ", " <<
        min_pxval[1] << ", " << min_pxval[2] << "}" << std::endl; std::cout <<
        "Per-channel MAX: {" << max_pxval[0] << ", " << max_pxval[1] << ", " <<
        max_pxval[2] << "}" << std::endl; pixel_t dynamic_channel_range =
        max_pxval - min_pxval; cv::divide(255, dynamic_channel_range,
        dynamic_channel_range); std::for_each(fflame_frame->data,
        fflame_frame->data + fflame_frame->rows * fflame_frame->cols,
                [min_pxval, dynamic_channel_range](pixel_t& pxval) {
                    pxval = (pxval - min_pxval).mul(dynamic_channel_range);
                });
        */
        cv::Mat_<pixel_t> keyframe_img(fflame_frame->rows, fflame_frame->cols,
                                       fflame_frame->data);

        // NOTE: we'll 'lose' the last frame generated in this manner?
        auto output_keyframe_fpath = output_dir;
        output_keyframe_fpath /= "baseframe_" + std::to_string(frame_idx);
        const std::string out_keyframe_filepath =
            output_keyframe_fpath.native();

        // write out the yml file...
        cv::FileStorage keyframe_ymlfile(out_keyframe_filepath + ".yml",
                                         cv::FileStorage::WRITE);
        keyframe_ymlfile << "baseframe_" + std::to_string(frame_idx)
                         << keyframe_img;

        cv::Mat logged_image;
        keyframe_img.convertTo(logged_image, CV_8UC3);
        //... and the image file
        cv::imwrite(out_keyframe_filepath + ".png", logged_image);
      }
    }
  }
}

template <template <class> class frame_t, typename pixel_t, typename data_t>
void generate_fractal_flames(const bfs::path output_dir,
                             const int render_height, const int render_width,
                             const int num_working_variants,
                             const int num_images,
                             const bool do_interpolation = true) {
  using flame_gen_t = fflame_generator<frame_t, data_t, pixel_t>;
  const int num_generator_threads = 4;
  auto bg_generator = std::unique_ptr<flame_gen_t>(
      new flame_gen_t(render_height, render_width, num_working_variants,
                      num_generator_threads));

  // pass in a queue to hold the finished flame frames
  // auto bg_framequeue = std::unique_ptr<EventQueue<frame_t<pixel_t>>>(new
  // EventQueue<frame_t<pixel_t>>(num_images));
  using bg_queue_t =
      typename flame_gen_t::template ts_queue_t<frame_t<pixel_t>>;
  bg_queue_t bg_framequeue;
  auto bg_framequeue_p = &bg_framequeue;
  bg_generator->register_framequeue(bg_framequeue_p);

  bg_generator->start_generation();
  run_frame_generation<bg_queue_t, frame_t, pixel_t>(
      bg_framequeue, output_dir, num_images, do_interpolation);
  bg_generator->stop_generation();
}

template <template <class> class frame_t, typename pixel_t, typename data_t>
void generate_fractal_flames(const bfs::path output_dir,
                             const int render_height, const int render_width,
                             std::vector<std::string> &&manual_variants,
                             const int num_images,
                             const bool do_interpolation = true) {
  using flame_gen_t = fflame_generator<frame_t, data_t, pixel_t>;
  const int num_generator_threads = 4;
  auto bg_generator = std::unique_ptr<flame_gen_t>(
      new flame_gen_t(render_height, render_width, std::move(manual_variants),
                      num_generator_threads));

  // pass in a queue to hold the finished flame frames
  // auto bg_framequeue = std::unique_ptr<EventQueue<frame_t<pixel_t>>>(new
  // EventQueue<frame_t<pixel_t>>(num_images));
  using bg_queue_t =
      typename flame_gen_t::template ts_queue_t<frame_t<pixel_t>>;
  bg_queue_t bg_framequeue;
  auto bg_framequeue_p = &bg_framequeue;
  bg_generator->register_framequeue(bg_framequeue_p);

  bg_generator->start_generation();
  run_frame_generation<bg_queue_t, frame_t, pixel_t>(
      bg_framequeue, output_dir, num_images, do_interpolation);
  bg_generator->stop_generation();
}

//------------------------------------------------------------------------------------------------------
template <typename pixel_t>
using frame_t = flame_frame<pixel_t>; // cv::Mat_<pixel_t>;

int main(int argc, char *argv[]) {
  std::vector<std::string> manual_variant_list;

  namespace bpo = boost::program_options;
  bpo::options_description bpo_desc("fractal flames options");
  bpo_desc.add_options()("help,h", "Print help message")(
      "output-path,o", bpo::value<std::string>(), "output flame-image path")(
      "num-variants,n", bpo::value<int>(),
      "number of working variants")("total-frames,t", bpo::value<int>(),
                                    "total number of keyframes to generate")(
      "variant-list,v",
      bpo::value<std::vector<std::string>>(&manual_variant_list),
      "list of variants to use")("do-interpolation,i", bpo::value<int>(),
                                 "flag for writing out interpolated frames "
                                 "between keyframes (1: yes, 0: no)");

  bpo::variables_map vm;
  bpo::store(bpo::command_line_parser(argc, argv).options(bpo_desc).run(), vm);
  bpo::notify(vm);

  std::string output_path;
  int num_working_variants, num_images, do_interpolation;

  if (vm.count("help")) {
    std::cout << bpo_desc << std::endl;
    return 0;
  }

  if (vm.count("output-path")) {
    output_path = vm["output-path"].as<std::string>();
  } else {
    std::cout << bpo_desc << std::endl;
    return 0;
  }

  if (vm.count("num-variants")) {
    num_working_variants = vm["num-variants"].as<int>();
  } else {
    std::cout << bpo_desc << std::endl;
    return 0;
  }

  if (vm.count("total-frames")) {
    num_images = vm["total-frames"].as<int>();
  } else {
    std::cout << bpo_desc << std::endl;
    return 0;
  }

  if (vm.count("do-interpolation")) {
    do_interpolation = vm["do-interpolation"].as<int>();
  } else {
    std::cout << bpo_desc << std::endl;
    return 0;
  }

  using data_t = double;
  using pixel_t = cv::Vec<data_t, 3>;
  const int render_height = 1024;
  const int render_width = 1024;

  const bool use_manual_variants = num_working_variants == 0;
  if (use_manual_variants) {
    // invoke via something like
    //./fflame_gen -o ff_frames/ff_v13 -n 10 -t 100 -i 0 -v fisheye -v eyefish
    //-v swirl -v this -v is -v a -v drill
    //
    //... where the -v are the manually provided variants, given in sequential
    //order next we just have to figure out how to pass these in...
    std::cout << "Given Variant (" << manual_variant_list.size()
              << ") List: " << std::endl;
    for (auto variant_id : manual_variant_list) {
      // search the list of variant names to make sure the provided ones are
      // valid
      std::string manual_var = variant_id;
      std::transform(manual_var.begin(), manual_var.end(), manual_var.begin(),
                     [](unsigned char c) { return std::tolower(c); });
      auto var_it = std::find(
          affine_fcns::variant_list<data_t>::variant_names.begin(),
          affine_fcns::variant_list<data_t>::variant_names.end(), manual_var);
      if (var_it != affine_fcns::variant_list<data_t>::variant_names.end()) {
        std::cout << variant_id << std::endl;
      } else {
        std::cout << "ERROR -- variant " << manual_var
                  << " is not valid. Valid variants list: " << std::endl;
        for (auto var_name : affine_fcns::variant_list<data_t>::variant_names) {
          std::cout << var_name << std::endl;
        }
        return 0;
      }
    }
  }

  // check the output directory path (if it doesn't exist, create it)
  bfs::path output_dir(output_path);
  if (!(bfs::exists(output_dir) && bfs::is_directory(output_dir))) {
    if (!boost::filesystem::create_directory(output_path)) {
      std::cout << "Invalid output file directory -- " << output_path
                << std::endl;
      throw std::invalid_argument("Invalid output file directory " +
                                  output_path);
    }
  }

  bool do_interpolation_flag = do_interpolation != 0;
  if (use_manual_variants) {
    generate_fractal_flames<frame_t, pixel_t, data_t>(
        output_dir, render_height, render_width, std::move(manual_variant_list),
        num_images, do_interpolation_flag);
  } else {
    generate_fractal_flames<frame_t, pixel_t, data_t>(
        output_dir, render_height, render_width, num_working_variants,
        num_images, do_interpolation_flag);
  }

  return 0;
}
