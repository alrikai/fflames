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

template <template <class> class frame_t, typename data_t, typename pixel_t = cv::Vec<data_t, 3>>
void render_fractal_flame(frame_t<pixel_t>* image, std::unique_ptr<std::vector<histogram_info<pixel_t>>> hist_data)
{
	const pixel_t background_pixel (0, 0, 0);
    frame_t<pixel_t> raw_image (fflame_constants::imheight, fflame_constants::imwidth, background_pixel);

    const int max_freq = (*std::max_element(hist_data->begin(), hist_data->end(), 
                [](histogram_info<pixel_t> lhs, histogram_info<pixel_t> rhs)
                { return lhs.frequency_count < rhs.frequency_count; })).frequency_count;
    //get the log of the maximum frequency count
    const data_t freq_max_log = std::log10(static_cast<data_t>(max_freq));
    
    //density estimation parameters
    const double min_est = 0.0;
    const double est_radius = 9.0;
    const double est_curve = 0.4;

    std::map<size_t, cv::Mat_<data_t>> density_est_kernels;
    cv::Mat_<uint64_t> image_density = cv::Mat_<uint64_t>::zeros(fflame_constants::imheight, fflame_constants::imwidth);

    //NOTE: we assume the histogram dimensions are multiples of the image dimensions
    const int rowpx_factor = fflame_constants::hist_height/fflame_constants::imheight;
    const int colpx_factor = fflame_constants::hist_width/fflame_constants::imwidth;

    int hist_data_idx = 0;
    //make the final output image -- NOTE: this is currently broken for any supersampling. 
    //we need to take both dimensions into account
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
                    auto h_data = hist_data->at(hist_data_idx++);
                    freq_avg += h_data.frequency_count;
                    color_avg += h_data.color;
                }
            }

            const auto freq_count = freq_avg;
            freq_avg /= rowpx_factor*colpx_factor;
            color_avg /= rowpx_factor*colpx_factor;

            data_t alpha = std::log10(freq_avg)/freq_max_log;
            raw_image.at(im_row, im_col) = 255 * color_avg * std::pow(alpha, fflame_constants::gamma_factor); 

            //next, we apply the density estimation filtering
            int kernel_width = std::lround(std::max(min_est, (est_radius / (std::pow(static_cast<double>(freq_count), est_curve)))));

            //dont do anything if we have a kernel size not-larger than 1 element
            if(kernel_width > 1)
            {
                //make sure we have an odd kernel size
                if(kernel_width % 2 == 0)
                    kernel_width++;

                //add the kernel matrix if it doesnt already exist
                auto kernel_it = density_est_kernels.find(kernel_width);
                if(kernel_it == density_est_kernels.end())
                {
                    //opencv default sigma: 0.3*((ksize-1)*0.5 - 1) + 0.8
                    cv::Mat_<data_t> kernel = cv::getGaussianKernel(kernel_width, -1);

                    cv::Mat_<data_t> kernel_2d = kernel*kernel.t();
                    //store a square 2D gaussian kernel
                    density_est_kernels.insert(std::make_pair(kernel_width, kernel_2d));
                }

                //associate the kernel to the current pixel
                image_density(im_row, im_col) = kernel_width;
            }
        }
    }

    //get rid of any NaNs
	for (int row = 0; row < raw_image.rows; ++row) {
        auto raw_rowp = raw_image.ptr(row);
		for (int col = 0; col < raw_image.cols; ++col, ++raw_rowp) {
            if(*raw_rowp != *raw_rowp) {
                *raw_rowp = background_pixel;
			} 
		}
	}
	
	/*
    cv::Mat_<cv::Vec<uint8_t, 3>> mask = cv::Mat(raw_image != raw_image);
    for (int r = 0; r < fflame_constants::imheight; ++r)
        for (int c = 0; c < fflame_constants::imwidth; ++c)
            for (int ch = 0; ch < 3; ++ch)
                if(mask(r,c)(ch))
                    raw_image(r,c)(ch) = 0;
    */

    //NOTE: have to apply the density estimation afterwards, since it requires pixels following the anchor pixel
    for (int im_row = 0; im_row < fflame_constants::imheight; ++im_row)
    {
        for (int im_col = 0; im_col < fflame_constants::imwidth; ++im_col)
        {
            const int kernel_width = image_density(im_row, im_col);
            const int kernel_hwidth = std::floor(kernel_width/2);
            auto kernel_it = density_est_kernels.find(kernel_width);
            if(kernel_it != density_est_kernels.end())
            {
                auto kernel = kernel_it->second;
                int kernel_ridx = 0;
                int kernel_cidx = 0;
                //clamps out of range indices at the borders
                for (int k_row = std::max(0, im_row-kernel_hwidth); k_row < std::min(fflame_constants::imheight-1, im_row+kernel_hwidth); ++k_row, ++kernel_ridx) {
                    for (int k_col = std::max(0, im_col-kernel_hwidth); k_col < std::min(fflame_constants::imwidth-1, im_col+kernel_hwidth); ++k_col, ++kernel_cidx) {
                        if (kernel_ridx < kernel.rows || kernel_cidx < kernel.cols) {
                            image->at(im_row, im_col) += raw_image.at(k_row, k_col) * kernel(kernel_ridx, kernel_cidx);
						} else {
                            std::cout << "NOTE: out of bounds on the kernel" << std::endl;
						}
					}
				}
            }
            else {
                image->at(im_row, im_col) = raw_image.at(im_row, im_col);
			}
        }
    }
}

/*
template <typename data_t, typename pixel_t>
void render_fractal_flame(cv::Mat_<pixel_t>& image)
{
    cv::Mat_<pixel_t> raw_image = cv::Mat_<data_t>::zeros(fflame_constants::imheight, fflame_constants::imwidth);
    auto& freq_hist = fflame_histdata<data_t, pixel_t>::fdata.freq_hist;
    auto& color_hist = fflame_histdata<data_t, pixel_t>::fdata.color_hist;

    const data_t freq_max_log = std::log10(*std::max_element(freq_hist.begin(), freq_hist.end()));
    
    //density estimation parameters
    const double min_est = 0.0;
    const double est_radius = 9.0;
    const double est_curve = 0.4;

    std::map<size_t, cv::Mat_<double>> density_est_kernels;
    cv::Mat_<uint64_t> image_density = cv::Mat_<uint64_t>::zeros(fflame_constants::imheight, fflame_constants::imwidth);

    //NOTE: we assume the histogram dimensions are multiples of the image dimensions
    const int rowpx_factor = fflame_constants::hist_height/fflame_constants::imheight;
    const int colpx_factor = fflame_constants::hist_width/fflame_constants::imwidth;
    //make the final output image
    for (int im_row = 0; im_row < fflame_constants::imheight; ++im_row) {
        int hist_row = im_row*rowpx_factor;     
        for (int im_col = 0; im_col < fflame_constants::imwidth; ++im_col) {
            data_t freq_avg = 0;
            pixel_t color_avg = 0;
            int hist_col = im_col*colpx_factor;    
            for (int supersample_row = 0; supersample_row < rowpx_factor; ++supersample_row) {
                for (int supersample_col = 0; supersample_col < colpx_factor; ++supersample_col) {
                    freq_avg += freq_hist(hist_row + supersample_row, hist_col + supersample_col);
                    color_avg += color_hist(hist_row + supersample_row, hist_col + supersample_col);
                }
            }

            const auto freq_count = freq_avg;
            freq_avg /= rowpx_factor*colpx_factor;
            color_avg /= rowpx_factor*colpx_factor;

            data_t alpha = std::log10(freq_avg)/freq_max_log;
            raw_image(im_row, im_col) = 255 * color_avg * std::pow(alpha, fflame_constants::gamma_factor); 

            //next, we apply the density estimation filtering
            int kernel_width = std::lround(std::max(min_est, (est_radius / (std::pow(static_cast<double>(freq_count), est_curve)))));

            //dont do anything if we have a kernel size not-larger than 1 element
            if(kernel_width > 1) {
                //make sure we have an odd kernel size
                if(kernel_width % 2 == 0) {
                    kernel_width++;
                }

                //add the kernel matrix if it doesnt already exist
                auto kernel_it = density_est_kernels.find(kernel_width);
                if(kernel_it == density_est_kernels.end()) {
                    //opencv default sigma: 0.3*((ksize-1)*0.5 - 1) + 0.8
                    cv::Mat_<double> kernel = cv::getGaussianKernel(kernel_width, -1);

                    cv::Mat_<double> kernel_2d = kernel*kernel.t();
                    //store a square 2D gaussian kernel
                    density_est_kernels.insert(std::make_pair(kernel_width, kernel_2d));
                    std::cout << "adding density estimation kernel @width " << kernel_width << " @ " << kernel_2d.rows << ", " << kernel_2d.cols << std::endl;
                }

                //associate the kernel to the current pixel
                image_density(im_row, im_col) = kernel_width;
            }
        }
    }

    //get rid of any NaNs
    cv::Mat_<cv::Vec<uint8_t, 3>> mask = cv::Mat(raw_image != raw_image);
    for (int r = 0; r < fflame_constants::imheight; ++r) {
        for (int c = 0; c < fflame_constants::imwidth; ++c) {
            for (int ch = 0; ch < 3; ++ch) {
                if(mask(r,c)(ch)) {
                    raw_image(r,c)(ch) = 0;
                }
            }
        }
    }

    std::cout << "density estimation: " << density_est_kernels.size() << " #kernels" << std::endl;

    //NOTE: have to apply the density estimation afterwards, since it requires pixels following the anchor pixel
    for (int im_row = 0; im_row < fflame_constants::imheight; ++im_row) {
        for (int im_col = 0; im_col < fflame_constants::imwidth; ++im_col) {
            const int kernel_width = image_density(im_row, im_col);
            const int kernel_hwidth = std::floor(kernel_width/2);
            auto kernel_it = density_est_kernels.find(kernel_width);
            if(kernel_it != density_est_kernels.end()) {
                auto kernel = kernel_it->second;
                int kernel_ridx = 0;
                int kernel_cidx = 0;
                //clamps out of range indices at the borders
                for (int k_row = std::max(0, im_row-kernel_hwidth); k_row < std::min(fflame_constants::imheight-1, im_row+kernel_hwidth); ++k_row, ++kernel_ridx) {
                    for (int k_col = std::max(0, im_col-kernel_hwidth); k_col < std::min(fflame_constants::imwidth-1, im_col+kernel_hwidth); ++k_col, ++kernel_cidx) {
                        image(im_row, im_col) += raw_image(k_row, k_col) * kernel(kernel_ridx, kernel_cidx);
                    }
                }
            }
            else {
                image(im_row, im_col) = raw_image(im_row, im_col);
            }
        }
    }

    static int counter = 0;
    const std::string raw_impath = "raw_image__" + std::to_string(counter++) + ".png";
    cv::imwrite(raw_impath, raw_image);
}
*/

#endif
