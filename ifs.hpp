#ifndef FF_IFS_HPP
#define FF_IFS_HPP

#include "ifs_types.hpp"

#include <cmath>
#include <cassert>

#include <string>
#include <vector>
#include <tuple>
#include <map>

#include <iostream>
#include <random>

namespace affine_fcns
{
    namespace detail 
    {
        template <typename data_t>
        void blend_colors(const std::vector<data_t>& color, flame_point<data_t>& point)
        {
            point.color[0] = (point.color[0] + color[0])/2.0;
            point.color[1] = (point.color[1] + color[1])/2.0;
            point.color[2] = (point.color[2] + color[2])/2.0;
        }
    }

    //just a linear function
    template <typename data_t>
    void V0 (flame_point<data_t>& point)
    {
        const static std::string name {"linear"};
        const static std::vector<data_t> color {0.25, 0.25, 1.0};
        detail::blend_colors(color, point);
    }

    //sinusoidal
    template <typename data_t>
    void V1 (flame_point<data_t>& point)
    {
        const static std::string name {"sinusoidal"};
        const static std::vector<data_t> color {1.0, 0.0, 0.0};

        point.x = std::sin(point.x);
        point.y = std::sin(point.y);
        detail::blend_colors(color, point);
    }

    //spherical
    template <typename data_t>
    void V2 (flame_point<data_t>& point)
    {
        const static std::string name {"spherical"};
        const static std::vector<data_t> color {0, 1.0, 0};        
        
        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        point.x /= r_factor;
        point.y /= r_factor;
        detail::blend_colors(color, point);
    }

    //swirl
    template <typename data_t>
    void V3 (flame_point<data_t>& point)
    {
        const static std::string name {"swirl"};
        const static std::vector<data_t> color {0.5, 0.1, 0.75};        
        
        auto theta = point.x * point.x + point.y * point.y;
        point.x = point.x * std::sin(theta) - point.y * std::cos(theta);
        point.y = point.x * std::cos(theta) + point.y * std::sin(theta);
        detail::blend_colors(color, point);
    }

    //horseshoe
    template <typename data_t>
    void V4 (flame_point<data_t>& point)
    {
        const static std::string name {"horseshoe"};
        const static std::vector<data_t> color {0, 0.75, 0.75};        

        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
     
        point.x = (point.x - point.y)*(point.x + point.y) / r_factor;
        point.y = 2 * point.x * point.y; 
        detail::blend_colors(color, point);
    }

    //polar
    template <typename data_t>
    void V5 (flame_point<data_t>& point)
    {   
        const static std::string name {"polar"};
        const static std::vector<data_t> color {0.5, 0.5, 0.5};        
    
        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);

        point.x = theta / fflame_constants::PI;
        point.y = r_factor - 1;

        detail::blend_colors(color, point);
    }

    //handkerchief
    template <typename data_t>
    void V6 (flame_point<data_t>& point)
    {   
        const static std::string name {"handkerchief"};
        const static std::vector<data_t> color {0.4, 0, 0.6};        
 
        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);
       
        point.x = r_factor * std::sin(theta + r_factor);
        point.y = r_factor * std::cos(theta - r_factor);

        detail::blend_colors(color, point);
    }


    //heart 
    template <typename data_t>
    void V7 (flame_point<data_t>& point)
    {   
        const static std::string name {"heart"};
        const static std::vector<data_t> color {0, 1.0, 0.3};        
 
        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);
       
        point.x = r_factor * std::sin(theta * r_factor);
        point.y = -r_factor * std::cos(theta * r_factor);       

        detail::blend_colors(color, point);
    }

    //disk
    template <typename data_t>
    void V8 (flame_point<data_t>& point)
    {   
        const static std::string name {"disk"};
        const static std::vector<data_t> color {0.9, 0.2, 0.2};        
 
        auto factor = fflame_constants::PI * std::sqrt(point.x * point.x + point.y * point.y);
        auto coeff = std::atan2(point.x, point.y) / fflame_constants::PI;

        point.x = coeff * std::sin(factor);
        point.y = coeff * std::cos(factor);       

        detail::blend_colors(color, point);
    }

    //spiral
    template <typename data_t>
    void V9 (flame_point<data_t>& point)
    {   
        const static std::string name {"spiral"};
        const static std::vector<data_t> color {0.35, 0.25, 0.7};        
     
        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);

        point.x = (std::cos(theta) + std::sin(r_factor)) / r_factor;
        point.y = (std::sin(theta) - std::cos(r_factor)) / r_factor;
   
        detail::blend_colors(color, point);
    }

    //hyperbolic
    template <typename data_t>
    void V10 (flame_point<data_t>& point)
    {   
        const static std::string name {"hyperbolic"};
        const static std::vector<data_t> color {0.1, 1.0, 0.1};        
 
        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);

        point.x = std::sin(theta) / r_factor;
        point.y = r_factor * std::cos(theta);
      
        detail::blend_colors(color, point);
    }

    //diamond
    template <typename data_t>
    void V11 (flame_point<data_t>& point)
    {   
        const static std::string name {"diamond"};
        const static std::vector<data_t> color {0.15, 0.15, 1.0};        

        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);

        point.x = std::sin(theta) * std::cos(r_factor);
        point.y = std::cos(theta) * std::sin(r_factor);

        detail::blend_colors(color, point);
    }

    //Ex
    template <typename data_t>
    void V12 (flame_point<data_t>& point)
    {   
        const static std::string name {"ex"};
        const static std::vector<data_t> color {0.2, 0.4, 0.8};        

        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);
        auto p0 = std::sin(theta + r_factor);
        auto p1 = std::cos(theta - r_factor);

        point.x = r_factor * (p0*p0*p0 + p1*p1*p1);
        point.y = r_factor * (p0*p0*p0 - p1*p1*p1);

        detail::blend_colors(color, point);
    }

    //julia
    template <typename data_t>
    void V13 (flame_point<data_t>& point)
    {   
        const static std::string name {"julia"};
        const static std::vector<data_t> color {0.8, 0.6, 0.2};        

        //no, the double sqrt is not a typo
        auto r_factor = std::sqrt(std::sqrt(point.x * point.x + point.y * point.y));
        auto theta = std::atan2(point.x, point.y);

        //omega is "a random variable that's either 0 or pi"
        static std::uniform_int_distribution<> omega_dist(0, 1); 
        double omega = fflame_constants::PI * omega_dist(fflame_randutil::get_engine());

        point.x = r_factor * std::cos(theta/2.0 + omega);
        point.y = r_factor * std::sin(theta/2.0 + omega);

        detail::blend_colors(color, point);
    }

    //bent
    template <typename data_t>
    void V14 (flame_point<data_t>& point)
    {   
        const static std::string name {"bent"};
        const static std::vector<data_t> color {0.25, 0.75, 0.25};        
       
        if(point.y < 0)
            point.y /= 2;
        if(point.x < 0)
            point.x *= 2;

        detail::blend_colors(color, point);
    }

    //waves
    template <typename data_t>
    void V15 (flame_point<data_t>& point)
    {   
        const static std::string name {"waves"};
        const static std::vector<data_t> color {0.5, 0.625, 0.125};        

        point.x += point.param_b * std::sin(point.y / (point.param_c * point.param_c));
        point.y += point.param_e * std::sin(point.x / (point.param_f * point.param_f)); 

        detail::blend_colors(color, point);
    }

    //fisheye
    template <typename data_t>
    void V16 (flame_point<data_t>& point)
    {   
        const static std::string name {"fisheye"};
        const static std::vector<data_t> color {0.1, 1.0, 0.4};        

        auto r_factor = 2.0 / (std::sqrt(point.x * point.x + point.y * point.y) + 1);
        
        //swap the x and y coordinates
        point.x = r_factor * point.y;
        point.y = r_factor * point.x;

        detail::blend_colors(color, point);
    }

    //popcorn
    template <typename data_t>
    void V17 (flame_point<data_t>& point)
    {   
        const static std::string name {"popcorn"};
        const static std::vector<data_t> color {0.0, 0.65, 0.9};        

        point.x += point.param_c * std::sin(std::tan(3 * point.y)); 
        point.y += point.param_f * std::sin(std::tan(3 * point.x)); 

        detail::blend_colors(color, point);
    }
   
    //exponential
    template <typename data_t>
    void V18 (flame_point<data_t>& point)
    {   
        const static std::string name {"exponential"};
        const static std::vector<data_t> color {0.5, 0.9, 0.0};        

        const double exp_factor = std::exp(point.x - 1);
        point.x = exp_factor * std::cos(fflame_constants::PI * point.y);
        point.y = exp_factor * std::sin(fflame_constants::PI * point.x);

        detail::blend_colors(color, point);
    }
    
    //power
    template <typename data_t>
    void V19 (flame_point<data_t>& point)
    {   
        const static std::string name {"power"};
        const static std::vector<data_t> color {0.0, 0.85, 0.45};        

        const auto theta = std::atan2(point.x, point.y);
        auto r_factor = std::pow(std::sqrt(point.x * point.x + point.y * point.y), std::sin(theta));
        
        point.x = r_factor * std::cos(theta);
        point.y = r_factor * std::sin(theta);

        detail::blend_colors(color, point);
    }
    
    
    //TODO still have 30 to go...

    template <typename data_t>
    struct invoker
    {
        void invoke(const size_t fcn_idx, flame_point<data_t>& pt) const
        {
            //apply the selected function's parameters to the input point
            assert(fcn_idx < fcn.size());
            auto flame_function = fcn[fcn_idx];

            //just for fun -- randomize the parameter lists too. Will want to change this to use the fast_rand 
            // -- apparently this takes ~15% of the TOTAL execution time (according to callgrind)
            static std::uniform_int_distribution<> fcn_dis(0, fcn.size()-1); 
            auto param_preit = affine_preparameters.find(fcn[fcn_dis(fflame_randutil::get_engine())]);
            //auto param_preit = affine_preparameters.find(flame_function);
            if(param_preit != affine_preparameters.end())
                pt.apply_affine_params(param_preit->second);

            flame_function (pt);

            //apply the post-processing affine transform parameters
            auto param_postit = affine_postparameters.find(fcn[fcn_dis(fflame_randutil::get_engine())]);
            //auto param_postit = affine_postparameters.find(flame_function);
            if(param_postit != affine_postparameters.end())
                pt.apply_affine_params(param_postit->second);        
        }

        //assign random values between [min_param, max_param) to all affine pre-process parameters
        void randomize_parameters(const data_t min_param, const data_t max_param)
        {
            //take the easy way out and assign uniform weighting to all functions
            const double probability_weights = 1.0 / fcn.size();
            fcn_probabilities.resize(fcn.size());
            std::fill(fcn_probabilities.begin(), fcn_probabilities.end(), probability_weights);

            std::uniform_real_distribution<> dis(min_param, max_param);
            for (size_t fcn_idx = 0; fcn_idx < fcn.size(); ++fcn_idx)
            {
                //make the pre-parameters
                auto param_preit = affine_preparameters.find(fcn[fcn_idx]);
                if(param_preit != affine_preparameters.end())
                {
                    flame_fcn_params<data_t> affine_params {
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine())
                    };
                    param_preit->second = std::move(affine_params);
                }

                //make the post-parameters
                auto param_postit = affine_postparameters.find(fcn[fcn_idx]);
                if(param_postit != affine_postparameters.end())
                {
                    flame_fcn_params<data_t> affine_params {
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine()),
                        dis(fflame_randutil::get_engine())
                    };
                    param_postit->second = std::move(affine_params);
                }
            }
        }

        typedef void (*flame_fcn)(flame_point<data_t>&);
        static std::vector<flame_fcn> fcn;
        static std::vector<data_t> fcn_probabilities;
        static std::map<flame_fcn, flame_fcn_params<data_t>> affine_preparameters;
        static std::map<flame_fcn, flame_fcn_params<data_t>> affine_postparameters;
    };

    template <typename data_t>
    std::vector<typename invoker<data_t>::flame_fcn> invoker<data_t>::fcn {
        {V0}, {V1}, {V2}, {V3}, {V4},
        {V5}, {V6}, {V7}, {V8}, {V9},
        {V10}, {V11}, {V12}, {V13}, {V14},    
        {V15}, {V16}, {V17}, {V18}, {V19}
    };

    //have some default parameter values for the functions, and allow user-modification of parameters 
    template <typename data_t>
    std::map<typename invoker<data_t>::flame_fcn, flame_fcn_params<data_t>> invoker<data_t>::affine_preparameters {
        {V0, flame_fcn_params<data_t>(-0.9825, -0.18, -0.6864, 0.0171, -0.2815, 0.0771)},
        {V1, flame_fcn_params<data_t>(-0.50123, 0.434, -0.0644, -0.4471, -0.0937, -0.8847)},
        {V2, flame_fcn_params<data_t>(0.1839, -0.1289, -0.1805, 0.0478, 0.18, 0.4097)},
        {V3, flame_fcn_params<data_t>(0.2, -0.25, -0.223, 0.244, 0, 1.23)},
        {V4, flame_fcn_params<data_t>(-0.15, -1.26, -0.2336, 0.2422, 0, 0.4455)}
    };
        
    template <typename data_t>
    std::map<typename invoker<data_t>::flame_fcn, flame_fcn_params<data_t>> invoker<data_t>::affine_postparameters  
    {};

    template <typename data_t>
    std::vector<data_t> invoker<data_t>::fcn_probabilities
    {};
}


#endif
