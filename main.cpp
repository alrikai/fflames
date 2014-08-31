#include "fractal_flame.hpp"

#include <opencv2/opencv.hpp>
#include <string>
#include <vector>
#include <memory>


double get_mat_pixel(cv::Mat_<cv::Vec<double,3>>& image, int row, int col, int channel)
{
    auto px = image(row, col);
    return px(channel);
}

//generates num_frames frames in interpolating between the lhs and rhs frames
template <typename pixel_t>
int interpolate_frames(const cv::Mat_<pixel_t>& lhs_img, const cv::Mat_<pixel_t>& rhs_img, const int num_frames, const std::string out_filebasepath, int frame_idx)
{

    std::string out_filepath = out_filebasepath + std::to_string(frame_idx++) + ".png";
    cv::imwrite(out_filepath, lhs_img);

    cv::Mat_<pixel_t> pixel_diffimg = (rhs_img - lhs_img) / static_cast<double>(num_frames);    
    //note: this is just a reference to lhs_img, not a clone
    cv::Mat_<pixel_t> working_frame = lhs_img;

    for (int i = 0; i < num_frames; ++i)
    {
        working_frame += pixel_diffimg;

        std::string out_filepath = out_filebasepath + std::to_string(frame_idx++) + ".png";
        cv::imwrite(out_filepath, working_frame);
    }

    //at this point, working_frame should be 1 step off from rhs_img
    return frame_idx;
}

int main(int argc, char* argv[])
{
    using data_t = double;
    using pixel_t = cv::Vec<data_t, 3>;  

    if(argc < 3)
    {
        std::cout << "Invalid Arguments -- specify output basepath for image and #frames" << std::endl;
        return 1;
    }

    std::vector<std::string> cmdline_args (argv + 1, argv + argc + !argc);
    std::string out_filebasepath = cmdline_args.at(0);
    const int num_images = std::stoi(cmdline_args.at(1));
    const uint8_t num_working_variants = 5;

    //for seeding the flame thread rng 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint8_t> working_variant_rng(0, num_working_variants-1);
    std::uniform_int_distribution<uint8_t> total_variant_rng(0, affine_fcns::variant_list<data_t>::variant_names.size() - 1);

    //make the starting variants (with parameters)
    //auto working_variants = generate_variants<data_t>(num_working_variants, variant_rng, gen);

    std::vector<std::shared_ptr<affine_fcns::variant<data_t>>> working_variants (num_working_variants);
    affine_fcns::variant_list<data_t> variant_maker;
    for (int i = 0; i < num_working_variants; ++i) 
    {
        auto selected_variant = affine_fcns::variant_list<data_t>::variant_names[total_variant_rng(gen)];
        working_variants.at(i) = std::shared_ptr<affine_fcns::variant<data_t>>(variant_maker.flame_maker.create_variant(selected_variant));
        // std::unique_ptr<affine_fcns::variant<data_t>>(variant_maker.flame_maker.create_variant(selected_variant));
    }

    affine_fcns::invoker<data_t> flamer (std::move(working_variants));
    flamer.randomize_parameters(-2, 2);

    cv::Mat_<pixel_t> prev_image;
    int output_index = 0;
    for (int i = 0; i < num_images; i++)
    {
        cv::Mat_<pixel_t> image;
        generate_fractal_flame<data_t, pixel_t>(flamer, image, 4);
        if(image.empty())
            std::cout << "invalid image" << std::endl;

        //get rid of any NaNs
        cv::Mat_<cv::Vec<uint8_t, 3>> mask = cv::Mat(image != image);
        for (int r = 0; r < image.rows; ++r)
            for (int c = 0; c < image.cols; ++c)
                for (int ch = 0; ch < 3; ++ch)
                    if(mask(r,c)(ch))
                        image(r,c)(ch) = 0;

        std::cout << "Writing frames [" << output_index;
        if(i > 0)
            output_index = interpolate_frames(prev_image, image, 30, out_filebasepath, output_index);
        std::cout << ", " << output_index << "]" << std::endl;

        auto selected_variant = affine_fcns::variant_list<data_t>::variant_names[total_variant_rng(gen)];
        //replace a random variant (that's not the linear variant)
        int mod_idx = working_variant_rng(gen);
        flamer.fcn.at(mod_idx).reset(variant_maker.flame_maker.create_variant(selected_variant)); 
        flamer.randomize_parameters(-2, 2);

        std::cout << "using variant " << selected_variant << std::endl;

        prev_image = image;
    }

    return 0;
}

/*
 * @TODO: figure out how best to create the animations. Do we generate a still frame using a subset of 
 * the variants, then tweak the variants/parameters/weights/probabilities, then generate another still
 * frame, then interpolate between the frames to have a smooth motion?
 * If so, we will need to be able to modify the following parameters:
 * - variation list
 * - variation affine parameter (pre and post transforms)
 * - variation weights
 * - variation probabilities
 *
 *   Since we want the movements between frames to be gradual, we will want to carry most of the same
 *   parameters through between the still frames. Hence, we can either have a long-running invoker object,
 *   or have some state object that we use to initialize the invoker each time. 
 *
 *   Also, we should move the random number generator state to a thread wrapper object, such that the 
 *   object has both the thread object and any per-thread state information (e.g. random number seed vals,
 *   in this case a uint64_t[2] array)
 */
