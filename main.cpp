#include "fractal_flame.hpp"

#include <opencv2/opencv.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <memory>
#include <stdexcept>

namespace bfs = boost::filesystem;

double get_mat_pixel(cv::Mat_<cv::Vec<double,3>>& image, int row, int col, int channel)
{
    auto px = image(row, col);
    return px(channel);
}

//generates num_frames frames in interpolating between the lhs and rhs frames
template <typename pixel_t>
int interpolate_frames(const cv::Mat_<pixel_t>& lhs_img, const cv::Mat_<pixel_t>& rhs_img, const int num_frames, const bfs::path out_filebasepath, int frame_idx)
{
  auto output_fpath = out_filebasepath;
  output_fpath /= std::to_string(frame_idx++) + ".png";
  const std::string out_filepath = output_fpath.native();
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
  }

  //at this point, working_frame should be 1 step off from rhs_img
  return frame_idx;
}

inline void print_ffhelp()
{
  std::cout << "program options are: " << std::endl;
	std::cout << "output flame-image path (output-path,o) " << std::endl;
  std::cout << "number of working variants (num-variants,n)" << std::endl;
  std::cout << "total number of frames to generate (total-frames,t)" << std::endl;
}

int main(int argc, char* argv[])
{
  namespace bpo = boost::program_options; 
  bpo::options_description bpo_desc("fractal flames options"); 
  bpo_desc.add_options() 
      ("help,h", "Print help message")
			("output-path,o", bpo::value<std::string>(), "output flame-image path")
			("num-variants,n", bpo::value<int>(), "number of working variants")
      ("total-frames,t", bpo::value<int>(), "total number of frames to generate");

	bpo::variables_map vm;
  bpo::store(bpo::command_line_parser(argc, argv).options(bpo_desc).run(), vm);
  bpo::notify(vm);

	std::string output_path;
	int num_working_variants, num_images;

	if(vm.count("help")){
    print_ffhelp();
		return 0;
	}

	if(vm.count("output-path")){
    output_path = vm["output-path"].as<std::string>();
	} else {
    print_ffhelp();
		return 0;
	}

	if(vm.count("num-variants")){
    num_working_variants = vm["num-variants"].as<int>();
	} else {
    print_ffhelp();
		return 0;
	}

	if(vm.count("total-frames")){
    num_images = vm["total-frames"].as<int>();
	} else {
    print_ffhelp();
		return 0;
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
        if(i > 0) {
            output_index = interpolate_frames(prev_image, image, 30, output_dir, output_index);
				}
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

