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

