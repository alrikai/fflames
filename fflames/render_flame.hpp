/* render_flame.hpp -- part of the fractal flames implementation
 *
 * Copyright (C) 2015 Alrik Firl
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef FF_RENDER_FLAME_HPP
#define FF_RENDER_FLAME_HPP

#include "ifs_types.hpp"

#include <map>
#include <opencv2/opencv.hpp>

namespace FFlames {

template <typename data_t> class fflame_renderer {
public:
  fflame_renderer(const int render_height, const int render_width,
                  const double min_est = 0.0, const double est_radius = 9.0,
                  const double est_curve = 0.4)
      : min_est(min_est), est_radius(est_radius), est_curve(est_curve),
        rowpx_factor(fflame_constants::imheight / render_height),
        colpx_factor(fflame_constants::imwidth / render_width) {}

  template <template <class> class frame_t,
            typename pixel_t = cv::Vec<data_t, 3>>
  void render(frame_t<pixel_t> *image,
              std::unique_ptr<std::vector<histogram_info<pixel_t>>> hist_data) {
    const pixel_t background_pixel(0, 0, 0);
    frame_t<pixel_t> raw_image(fflame_constants::imheight,
                               fflame_constants::imwidth, background_pixel);
    const int max_freq =
        (*std::max_element(
             hist_data->begin(), hist_data->end(),
             [](histogram_info<pixel_t> lhs, histogram_info<pixel_t> rhs) {
               return lhs.frequency_count < rhs.frequency_count;
             }))
            .frequency_count;
    // get the log of the maximum frequency count
    const data_t freq_max_log = std::log10(static_cast<data_t>(max_freq));
    compute_rawimage_density<frame_t, pixel_t>(raw_image, freq_max_log,
                                               std::move(hist_data));

    //--------------------------------------------------------------------------------------------
    // get rid of any NaNs
    for (int row = 0; row < raw_image.rows; ++row) {
      auto raw_rowp = raw_image.ptr(row);
      for (int col = 0; col < raw_image.cols; ++col, ++raw_rowp) {
        auto color_px = *raw_rowp;
        for (int color_idx = 0; color_idx < 3; color_idx++) {
          if (std::isnan(color_px[color_idx]) ||
              std::isinf(color_px[color_idx])) {
            std::cout << "NOTE: got NaN/Inf in rendered image" << std::endl;
          }
        }
      }
    }
    //--------------------------------------------------------------------------------------------

    render_frame_density<frame_t, pixel_t>(std::move(raw_image), image);
  }

private:
  template <template <class> class frame_t,
            typename pixel_t = cv::Vec<data_t, 3>>
  void compute_rawimage_density(
      frame_t<pixel_t> &raw_image, const data_t freq_max_log,
      std::unique_ptr<std::vector<histogram_info<pixel_t>>> hist_data);

  template <template <class> class frame_t,
            typename pixel_t = cv::Vec<data_t, 3>>
  void render_frame_density(frame_t<pixel_t> &&raw_image,
                            frame_t<pixel_t> *image);

  std::map<size_t, cv::Mat_<data_t>> density_est_kernels;
  cv::Mat_<uint64_t> image_density;

  // density estimation parameters
  const double min_est;
  const double est_radius;
  const double est_curve;

  // NOTE: we assume the histogram dimensions are multiples of the image
  // dimensions
  const int rowpx_factor;
  const int colpx_factor;
};

template <typename data_t>
template <template <class> class frame_t, typename pixel_t>
void fflame_renderer<data_t>::compute_rawimage_density(
    frame_t<pixel_t> &raw_image, const data_t freq_max_log,
    std::unique_ptr<std::vector<histogram_info<pixel_t>>> hist_data) {
  // reset the image density map for this frame
  image_density = cv::Mat_<uint64_t>::zeros(fflame_constants::imheight,
                                            fflame_constants::imwidth);

  int hist_data_idx = 0;
  // make the final output image -- NOTE: this is currently broken for any
  // supersampling. we need to take both dimensions into account
  for (int im_row = 0; im_row < fflame_constants::imheight; ++im_row) {
    // int hist_row = im_row*rowpx_factor;
    for (int im_col = 0; im_col < fflame_constants::imwidth; ++im_col) {
      data_t freq_count = 0;
      pixel_t color_avg = 0;
      // int hist_col = im_col*colpx_factor;
      for (int supersample_row = 0; supersample_row < rowpx_factor;
           ++supersample_row) {
        for (int supersample_col = 0; supersample_col < colpx_factor;
             ++supersample_col) {
          auto h_data = hist_data->at(hist_data_idx++);
          freq_count += h_data.frequency_count;
          color_avg += h_data.color;
        }
      }

      //--------------------------------------------------------------------------------------------
      // check for conditions that cause numerical instability
      if (std::abs(freq_count) < 1.00001f ||
          (color_avg[0] + color_avg[1] + color_avg[2] == 0)) {
        raw_image.at(im_row, im_col) = color_avg;
        continue;
      }
      //--------------------------------------------------------------------------------------------

      /*
      const auto freq_count = freq_avg;
      freq_avg /= rowpx_factor*colpx_factor;
      color_avg /= rowpx_factor*colpx_factor;
      data_t alpha = std::log10(freq_avg)/freq_max_log;

      //can use this to adjust the brightness
      static constexpr data_t color_pixelfactor = 200;
      auto color_px = color_pixelfactor * color_avg * std::pow(alpha,
      fflame_constants::gamma_factor);
      */

      const auto freq_avg = freq_count / rowpx_factor * colpx_factor;
      color_avg /= rowpx_factor * colpx_factor;
      const data_t alpha = std::log10(freq_avg) / freq_avg;

      // can use this to adjust the brightness
      const data_t color_pixelfactor = 128 * alpha;
      auto px_gammacorrected =
          std::pow(freq_avg, fflame_constants::gamma_factor);
      auto color_px = color_pixelfactor * color_avg * px_gammacorrected;

      //--------------------------------------------------------------------------------------------
      for (int color_idx = 0; color_idx < 3; color_idx++) {
        if (std::isnan(color_px[color_idx]) ||
            std::isinf(color_px[color_idx])) {
          std::cout << "NOTE: got NaN/Inf in rendered image" << std::endl;
        }
      }
      //--------------------------------------------------------------------------------------------
      raw_image.at(im_row, im_col) = color_px;

      // next, we apply the density estimation filtering
      int kernel_width = std::lround(std::max(
          min_est, (est_radius /
                    (std::pow(static_cast<double>(freq_count), est_curve)))));

      // dont do anything if we have a kernel size not-larger than 1 element
      if (kernel_width > 1) {
        // make sure we have an odd kernel size
        if (kernel_width % 2 == 0) {
          kernel_width++;
        }

        // add the kernel matrix if it doesnt already exist
        auto kernel_it = density_est_kernels.find(kernel_width);
        if (kernel_it == density_est_kernels.end()) {
          // opencv default sigma: 0.3*((ksize-1)*0.5 - 1) + 0.8
          cv::Mat_<data_t> kernel = cv::getGaussianKernel(kernel_width, -1);
          cv::Mat_<data_t> kernel_2d = kernel * kernel.t();
          // store a square 2D gaussian kernel
          density_est_kernels.insert(std::make_pair(kernel_width, kernel_2d));
        }

        // associate the kernel to the current pixel
        image_density(im_row, im_col) = kernel_width;
      }
    }
  }
}

template <typename data_t>
template <template <class> class frame_t, typename pixel_t>
void fflame_renderer<data_t>::render_frame_density(frame_t<pixel_t> &&raw_image,
                                                   frame_t<pixel_t> *image) {

  auto max_possible_value =
      std::numeric_limits<typename pixel_t::value_type>::max();
  auto min_possible_value =
      std::numeric_limits<typename pixel_t::value_type>::min();
  pixel_t min_pixel(max_possible_value, max_possible_value, max_possible_value);
  pixel_t max_pixel(min_possible_value, min_possible_value, min_possible_value);

  // NOTE: have to apply the density estimation afterwards, since it requires
  // pixels following the anchor pixel
  for (int im_row = 0; im_row < fflame_constants::imheight; ++im_row) {
    for (int im_col = 0; im_col < fflame_constants::imwidth; ++im_col) {
      const int kernel_width = image_density(im_row, im_col);
      const int kernel_hwidth = std::floor(kernel_width / 2);
      auto kernel_it = density_est_kernels.find(kernel_width);
      if (kernel_it != density_est_kernels.end()) {
        auto kernel = kernel_it->second;
        int kernel_ridx = 0;
        int kernel_cidx = 0;
        // clamps out of range indices at the borders
        for (int k_row = std::max(0, im_row - kernel_hwidth);
             k_row <
             std::min(fflame_constants::imheight - 1, im_row + kernel_hwidth);
             ++k_row, ++kernel_ridx) {
          for (int k_col = std::max(0, im_col - kernel_hwidth);
               k_col <
               std::min(fflame_constants::imwidth - 1, im_col + kernel_hwidth);
               ++k_col, ++kernel_cidx) {
            if (kernel_ridx < kernel.rows || kernel_cidx < kernel.cols) {
              image->at(im_row, im_col) +=
                  raw_image.at(k_row, k_col) * kernel(kernel_ridx, kernel_cidx);
            } else {
              std::cout << "NOTE: out of bounds on the kernel" << std::endl;
            }
          }
        }
      } else {
        image->at(im_row, im_col) = raw_image.at(im_row, im_col);
      }

      for (int k = 0; k < 3; k++) {
        if (image->at(im_row, im_col)[k] > max_pixel[k]) {
          max_pixel[k] = image->at(im_row, im_col)[k];
        }
        if (image->at(im_row, im_col)[k] < min_pixel[k]) {
          min_pixel[k] = image->at(im_row, im_col)[k];
        }
      }
    }
  }

  // Q: at this point, the data will no longer be within the [0, 256) range.
  // Should we re-normalize here?
  //--> or should we let the caller do so (i.e. if they want 16bpp images?)...
  // might as well re-normalize here

  /*
  //find the max pixel value... the question is, how to do this w/ RGB? --> per
  channel basis, for now const pixel_t max_norm_pixel_value (255 / (max_pixel[0]
  - min_pixel[0]), 255 / (max_pixel[1] - min_pixel[1]), 255 / (max_pixel[2] -
  min_pixel[2])); for (int im_row = 0; im_row < fflame_constants::imheight;
  ++im_row)
  {
      for (int im_col = 0; im_col < fflame_constants::imwidth; ++im_col)
      {
          image->at(im_row, im_col) = (image->at(im_row, im_col) -
  min_pixel).mul(max_norm_pixel_value);
      }
  }
  */
}

}

#endif
