/* fflame_generator.hpp -- part of the fractal flames implementation
 *
 * Copyright (C) 2015 Alrik Firl
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef FF_FFLAME_GENERATOR_HPP
#define FF_FFLAME_GENERATOR_HPP

#include "fractal_flame.hpp"
#include "ifs.hpp"
#include "ifs_types.hpp"
#include "render_flame.hpp"
#include "util/EventQueue.hpp"
#include "util/ff_utils.hpp"

#include <opencv2/opencv.hpp>
//#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>

#include <atomic>
#include <condition_variable>
#include <limits>
#include <memory>
#include <mutex>
#include <random>
#include <string>
#include <thread>

namespace FFlames {

template <template <class> class frame_t, typename data_t, typename pixel_t>
class fflame_generator {
public:
  // NOTE: can (try to) get a reproducable sequence by not initializing
  // flame_gen(flame_rd())?
  fflame_generator(const int imheight, const int imwidth,
                   const int num_variants,
                   const int num_workers = std::thread::hardware_concurrency())
      : num_workers(num_workers), render_imheight(imheight),
        render_imwidth(imwidth), num_working_variants(num_variants),
        fflame_histoqueue(100, 30), flame_prebarrier(num_workers),
        flame_postbarrier(num_workers), fflame_th(nullptr),
        fflame_histdata(nullptr), flame_gen(flame_rd()),
        total_variant_rng(
            0, affine_fcns::variant_list<data_t>::variant_names.size() - 1) {
    // NOTE: this has to be shared between all the worker threads
    fflame_histdata = std::unique_ptr<fflame_data<data_t, pixel_t>>(
        new fflame_data<data_t, pixel_t>());

    initialize_variants();
    fflame_state.store(false);
  }

  fflame_generator(const int imheight, const int imwidth,
                   std::vector<std::string> &&manual_variants,
                   const int num_workers = std::thread::hardware_concurrency())
      : num_workers(num_workers), render_imheight(imheight),
        render_imwidth(imwidth), num_working_variants(manual_variants.size()),
        fflame_histoqueue(100, 30), flame_prebarrier(num_workers),
        flame_postbarrier(num_workers), fflame_th(nullptr),
        fflame_histdata(nullptr), flame_gen(flame_rd()),
        total_variant_rng(
            0, affine_fcns::variant_list<data_t>::variant_names.size() - 1) {
    // NOTE: this has to be shared between all the worker threads
    fflame_histdata = std::unique_ptr<fflame_data<data_t, pixel_t>>(
        new fflame_data<data_t, pixel_t>());

    initialize_variants(std::move(manual_variants));
    fflame_state.store(false);
  }

  ~fflame_generator() {
    if (fflame_state.load()) {
      stop_generation();
    }
  }

#if 1
  template <typename T> using ts_queue_t = EventQueue<T>;
#else
  template <typename T>
  using ts_queue_t =
      boost::lockfree::spsc_queue<T, boost::lockfree::capacity<1024>>;
#endif

  void start_generation() {
    // if already started, don't try to start again
    if (fflame_state.load()) {
      return;
    }

    fflame_state.store(true);
    fflame_th = std::unique_ptr<std::thread>(
        new std::thread(&fflame_generator::start_fflame_generation, this));
  }

  void stop_generation() {
    if (fflame_state.load()) {
      fflame_state.store(false);
      for (size_t th_idx = 0; th_idx < fflame_workers.size(); ++th_idx) {
        fflame_workers.at(th_idx)->finish_flame();
      }
      fflame_th->join();
    }
  }

  void register_framequeue(ts_queue_t<frame_t<pixel_t>> *queue) {
    fflame_imagequeue = queue;
  }

private:
  void initialize_variants(std::vector<std::string> &&manual_variants);
  void initialize_variants();
  void start_fflame_generation();
  void generate_fflame(fflame_util::fast_rand &rand_gen);
  void render_fflame();

  // number of threads used for the generation (not counting rendering)
  int num_workers;
  int render_imheight;
  int render_imwidth;

  // the number of variants to have active
  uint8_t num_working_variants;
  std::thread::id worker_overlord_id;

  // pause generation if above the max, resume if paused and below the min
  static constexpr size_t max_image_thresh = 25;
  static constexpr size_t min_image_thresh = 5;

  // controls the starting and stopping of the worker threads
  std::mutex gen_mtx;
  std::condition_variable gen_cv;

  // for passing results asynchronously between the generate & render steps
  // EventQueue<cv::Mat_<histogram_info<pixel_t>>> fflame_histoqueue;

  ts_queue_t<std::vector<histogram_info<pixel_t>>> fflame_histoqueue;
  // for holding the final images. provided by the caller
  ts_queue_t<frame_t<pixel_t>> *fflame_imagequeue;

  fflame_util::barrier flame_prebarrier;
  fflame_util::barrier flame_postbarrier;

  // the rendering/controller flame thread
  std::unique_ptr<std::thread> fflame_th;
  // the fractal flame generation threads
  std::vector<std::unique_ptr<flame_thread>> fflame_workers;

  // flag for starting/stopping the whole fflame pipeline
  std::atomic<bool> fflame_state;

  // holds the list of the current variation functions to use
  std::unique_ptr<affine_fcns::invoker<data_t>> flamer;
  std::unique_ptr<fflame_data<data_t, pixel_t>> fflame_histdata;

  // the various RNGs needed for the generation
  std::random_device flame_rd;
  std::mt19937 flame_gen;
  // std::uniform_int_distribution<uint8_t> working_variant_rng;
  std::uniform_int_distribution<> total_variant_rng;
};

template <template <class> class frame_t, typename data_t, typename pixel_t>
void fflame_generator<frame_t, data_t, pixel_t>::initialize_variants(
    std::vector<std::string> &&manual_variants) {
  flamer = std::unique_ptr<affine_fcns::invoker<data_t>>(
      new affine_fcns::invoker<data_t>(std::move(manual_variants)));
  flamer->randomize_parameters(-2, 2);
}

template <template <class> class frame_t, typename data_t, typename pixel_t>
void fflame_generator<frame_t, data_t, pixel_t>::initialize_variants() {
  flamer = std::unique_ptr<affine_fcns::invoker<data_t>>(
      new affine_fcns::invoker<data_t>(num_working_variants));
  // NOTE: we always want to have the linear variant (variant #0)
  int variant_idx = 0;
  auto linear_variant_id =
      affine_fcns::variant_list<data_t>::variant_names[variant_idx];
  flamer->set_variant(variant_idx, linear_variant_id);
  for (int i = variant_idx + 1; i < num_working_variants; ++i) {
    std::string selected_variant;
    // keep re-rolling the variant until we get a valid one
    do {
      variant_idx = total_variant_rng(flame_gen);
      selected_variant =
          affine_fcns::variant_list<data_t>::variant_names[variant_idx];
    } while (!flamer->check_variant_valid(selected_variant));

    // can assume by here that the selected variant is a valid one
    flamer->set_variant(i, selected_variant);
  }
  flamer->randomize_parameters(-2, 2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <template <class> class frame_t, typename data_t, typename pixel_t>
void fflame_generator<frame_t, data_t, pixel_t>::start_fflame_generation() {
  fflame_state.store(true);

  std::uniform_int_distribution<uint64_t> flame_dist(
      0, std::numeric_limits<uint64_t>::max());

  // spawn the worker threads
  for (int th_idx = 0; th_idx < num_workers; ++th_idx) {
    // seeds the flame thread's fast random number generator
    fflame_workers.emplace_back(
        new flame_thread(flame_dist(flame_gen), flame_dist(flame_gen)));
    // this step is the one that actually spawns the thread
    bool flame_run;
    std::thread::id tid;
    std::tie(flame_run, tid) = fflame_workers.at(th_idx)->do_flame(
        &fflame_generator::generate_fflame, this);

    // this would be sort of a big deal. should probably throw something here
    if (!flame_run) {
      std::cout << "NOTE: flame task failed" << std::endl;
    }

    if (th_idx == 0) {
      worker_overlord_id = tid;
    }
  }

  // prepare the image and push the image to the shared queue
  render_fflame();
}

template <typename pixel_t>
histogram_info<pixel_t>
print_mat(const cv::Mat_<histogram_info<pixel_t>> &hdata, int row, int col) {
  return hdata(row, col);
}

// makes the fractal flame histogram using the chaos game. Is invoked by N
// threads
template <template <class> class frame_t, typename data_t, typename pixel_t>
void fflame_generator<frame_t, data_t, pixel_t>::generate_fflame(
    fflame_util::fast_rand &rand_gen) {
  while (fflame_state.load()) {
    auto start_iter_time = std::chrono::high_resolution_clock::now();

    // 0. check if the generation should pause (i.e. too many frames in queue)
    // --> leave this for after you get it working
    // 1. start the next flame generation
    run_fflame<data_t, pixel_t>(flamer.get(),
                                fflame_constants::num_samples / num_workers,
                                fflame_histdata.get(), rand_gen);

    auto end_iter_time = std::chrono::high_resolution_clock::now();
    double time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                              end_iter_time - start_iter_time)
                              .count();
    std::cout << "thread " << std::this_thread::get_id()
              << " generation took -- " << time_elapsed << " ms" << std::endl;

    // 2. wait for all the threads to finish the previous round
    flame_prebarrier.wait();

    // merge the N thread results using the overlord thread
    if (std::this_thread::get_id() == worker_overlord_id) {
      static std::vector<std::string> mutated_variant_ids(num_working_variants,
                                                          "default");

      // somewhat unfortunate, but need to dynamically allocate to avoid scoping
      // problems
      auto hist_info = std::unique_ptr<std::vector<histogram_info<pixel_t>>>(
          new std::vector<histogram_info<pixel_t>>(fflame_constants::imheight *
                                                   fflame_constants::imwidth));
      fflame_histdata->get_and_reset(*hist_info);

      // 3.5 pass the finished histogram to the shared-buffer for rendering
      fflame_histoqueue.push(std::move(hist_info));

      // 4. mutate the variants
      std::string selected_variant;
      // keep re-rolling the variant until we get a valid one
      do {
        selected_variant =
            affine_fcns::variant_list<data_t>::variant_names[total_variant_rng(
                flame_gen)];
      } while (!flamer->check_variant_valid(selected_variant));

      // replace a random variant (that's not the linear variant)
      int mod_idx =
          (num_working_variants > 1)
              ? total_variant_rng(flame_gen) % (num_working_variants - 1) + 1
              : 0;
      // can assume by here that the selected variant is a valid one
      flamer->set_variant(mod_idx, selected_variant);
      flamer->randomize_parameters(-2, 2);

      flamer->print_variant_list();

      // there should be no threads waiting at the pre-barrier at this point;
      // the other (non-overlord) threads should be waiting at the post-barrier
      flame_prebarrier.reset(num_workers);
    }
  }
}

// generates an image based on the fflame histogram
template <template <class> class frame_t, typename data_t, typename pixel_t>
void fflame_generator<frame_t, data_t, pixel_t>::render_fflame() {
  bool got_histdata = false;
  // int raw_counter = 0;

  fflame_renderer<data_t> flame_image_renderer(render_imheight, render_imwidth);

  double total_render_time = 0;
  while (fflame_state.load()) {
    // 1. get the histogram
    auto hist_info = fflame_histoqueue.pop(got_histdata);
    if (got_histdata && hist_info) {
      auto start_render_time = std::chrono::high_resolution_clock::now();

      std::unique_ptr<frame_t<pixel_t>> image =
          std::unique_ptr<frame_t<pixel_t>>(
              new frame_t<pixel_t>(render_imheight, render_imwidth));
      std::fill(image->data, image->data + image->rows * image->cols, 0);

      // 2. call the rendering routine, get resultant image
      flame_image_renderer.template render<frame_t, pixel_t>(
          image.get(), std::move(hist_info));

      /*
                  //filter out the flames that are too sparse
                  int num_nonzero = 0;
                  const double image_threshold = 0.1 * imwidth * imheight;
                  const double px_threshold = 1.0;
                  for (int r = 0; r < imheight; ++r)
                  {
                      for (int c = 0; c < imwidth; ++c)
                      {
                          size_t px_sum = cv::sum(outfile_image.at<pixel_t>(r,
         c))[0]; if(px_sum > px_threshold) { num_nonzero++;
                          }
                      }
                  }
                  if(num_nonzero < image_threshold)
                  {
                      std::cout << "NOTE: flame too dark" << std::endl;
                      continue;
                  }
                  //mostly for debugging -- save the images to disk
                  const std::string raw_impath = "FFLAME_baseimage_" +
         std::to_string(raw_counter++) + ".png"; cv::imwrite(raw_impath, image);
      */

      fflame_imagequeue->push(std::move(image));

      auto end_render_time = std::chrono::high_resolution_clock::now();
      total_render_time +=
          std::chrono::duration_cast<std::chrono::milliseconds>(
              end_render_time - start_render_time)
              .count();
    }
  }

  std::cout << "rendering took -- " << total_render_time << " ms" << std::endl;
}

} // namespace FFlames
#endif
