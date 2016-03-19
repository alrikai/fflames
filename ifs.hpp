/* ifs.hpp -- part of the fractal flames implementation 
 *
 * Copyright (C) 2015 Alrik Firl 
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */


#ifndef FF_IFS_HPP
#define FF_IFS_HPP

#include "util/factory.hpp"
#include "ifs_types.hpp"

#include <cmath>
#include <cassert>

#include <string>
#include <vector>
#include <tuple>
#include <map>

#include <iostream>
#include <random>
#include <numeric>

namespace affine_fcns
{
    namespace detail 
    {
        template <typename data_t>
        inline void blend_colors(const std::vector<data_t>& color, flame_point<data_t>& point)
        {
            point.color[0] = (point.color[0] + color[0])/2.0;
            point.color[1] = (point.color[1] + color[1])/2.0;
            point.color[2] = (point.color[2] + color[2])/2.0;
        }

        constexpr double eps = 1e-10;
    }

    template <typename data_t>
    struct variant
    {
        variant(const data_t weight, const std::vector<data_t>& color)
            : weight(weight), color(color)
        {}

        virtual void apply_variant(flame_point<data_t>& point) = 0;

        double weight;
        std::vector<data_t> color;
    };
    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V0 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;
         
        explicit V0(const data_t weight, const std::vector<data_t>& color = std::vector<data_t>{0.25, 0.25, 1.0})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            point.x = weight * std::sin(point.x);
            point.y = weight * std::sin(point.y);
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V0<data_t>::name {"linear"};
    //////////////////////////////////////////////////////////////////

    //sinusoidal
    template <typename data_t>
    struct V1 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;
         
        explicit V1(const data_t weight, const std::vector<data_t>& color = std::vector<data_t>{1.0, 0.0, 0.0})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            point.x = weight * std::sin(point.x);
            point.y = weight * std::sin(point.y);
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V1<data_t>::name {"sinusoidal"};
    //////////////////////////////////////////////////////////////////

    //spherical
    template <typename data_t>
    struct V2 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V2(const data_t weight, const std::vector<data_t>& color = std::vector<data_t>{0.0, 1.0, 0.0})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = 1.0 / (point.x * point.x + point.y * point.y + detail::eps);
            point.x = weight * point.x * r_factor;
            point.y = weight * point.y * r_factor;
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V2<data_t>::name {"spherical"};
    //////////////////////////////////////////////////////////////////
    
    //swirl
    template <typename data_t>
    struct V3 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V3(const data_t weight, const std::vector<data_t>& color = std::vector<data_t>{0.5, 0.1, 0.75})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto theta = point.x * point.x + point.y * point.y;
            point.x = weight * point.x * std::sin(theta) - point.y * std::cos(theta);
            point.y = weight * point.x * std::cos(theta) + point.y * std::sin(theta);
            detail::blend_colors(color, point);        
        }
    };    
    template <typename data_t>
    const std::string V3<data_t>::name {"swirl"};
    //////////////////////////////////////////////////////////////////
   
    //horseshoe
    template <typename data_t>
    struct V4 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V4(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0, 0.75, 0.75})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = 1.0 / (std::sqrt(point.x * point.x + point.y * point.y) + detail::eps);
         
            point.x = weight * (point.x - point.y)*(point.x + point.y) * r_factor;
            point.y = weight * 2 * point.x * point.y * r_factor; 
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V4<data_t>::name {"horseshoe"};
    //////////////////////////////////////////////////////////////////

    //polar
    template <typename data_t>
    struct V5 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V5(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.5, 0.5, 0.5})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
            auto theta = std::atan2(point.x, point.y);

            point.x = weight * theta / fflame_constants::PI;
            point.y = weight * r_factor - 1.0;

            detail::blend_colors(color, point);       
        }
    };
    template <typename data_t>
    const std::string V5<data_t>::name {"polar"};
    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V6 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V6(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.4, 0.0, 0.6})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = weight * std::sqrt(point.x * point.x + point.y * point.y);
            auto theta = std::atan2(point.x, point.y);
             
            point.x = r_factor * std::sin(theta + r_factor);
            point.y = r_factor * std::cos(theta - r_factor);
      
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V6<data_t>::name {"handkerchief"};
    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V7 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V7(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0, 1.0, 0.3})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = weight * std::sqrt(point.x * point.x + point.y * point.y);
            auto theta = std::atan2(point.x, point.y);
           
            point.x = r_factor * std::sin(theta * r_factor);
            point.y = -r_factor * std::cos(theta * r_factor);       

            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V7<data_t>::name {"heart"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V8 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V8(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.9, 0.2, 0.2})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto factor = fflame_constants::PI * std::sqrt(point.x * point.x + point.y * point.y);
            auto coeff = weight * std::atan2(point.x * fflame_constants::PI, point.y * fflame_constants::PI) / fflame_constants::PI;
      
            point.x = coeff * std::sin(factor);
            point.y = coeff * std::cos(factor);       
      
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V8<data_t>::name {"disk"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V9 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V9(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.35, 0.25, 0.75})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
            auto theta = std::atan2(point.x, point.y);
      
            point.x = weight * (std::cos(theta) + std::sin(r_factor)) / (r_factor + detail::eps);
            point.y = weight * (std::sin(theta) - std::cos(r_factor)) / (r_factor + detail::eps);
         
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V9<data_t>::name {"spiral"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V10 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V10(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0, 0, 0})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
            auto theta = std::atan2(point.x, point.y);
 
            point.x = weight * std::sin(theta) / r_factor;
            point.y = weight * r_factor * std::cos(theta);
        
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V10<data_t>::name {"hyperbolic"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V11 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V11(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.15, 0.25, 1.0})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
            auto theta = std::atan2(point.x, point.y);
  
            point.x = weight * std::sin(theta) * std::cos(r_factor);
            point.y = weight * std::cos(theta) * std::sin(r_factor);
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V11<data_t>::name {"diamond"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V12 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V12(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.2, 0.4, 0.8})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
            auto theta = std::atan2(point.x, point.y);
            auto p0 = std::sin(theta + r_factor);
            auto p1 = std::cos(theta - r_factor);
  
            point.x = weight * r_factor * (p0*p0*p0 + p1*p1*p1);
            point.y = weight * r_factor * (p0*p0*p0 - p1*p1*p1);
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V12<data_t>::name {"ex"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V13 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V13(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.8, 0.4, 0.2})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            //no, the double sqrt is not a typo
            auto r_factor = weight * std::pow(std::sqrt(std::sqrt(point.x * point.x + point.y * point.y)), 0.25);
            auto theta = std::atan2(point.x, point.y) / 2.0;
  
            //omega is "a random variable that's either 0 or pi"
            static std::uniform_int_distribution<> omega_dist(0, 1); 
            theta += fflame_constants::PI * omega_dist(fflame_util::get_engine());
  
            point.x = r_factor * std::cos(theta);
            point.y = r_factor * std::sin(theta);
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V13<data_t>::name {"julia"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V14 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V14(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.25, 0.75, 0.25})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            if(point.y < 0)
                point.y /= 2;
            if(point.x < 0)
                point.x *= 2;
  
            point.x *= weight;
            point.y *= weight;
  
            detail::blend_colors(color, point); 
        }
    };
    template <typename data_t>
    const std::string V14<data_t>::name {"bent"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V15 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V15(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.625, 0.625, 0.125})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            point.x += point.param_b * std::sin(point.y / (point.param_c * point.param_c + detail::eps));
            point.y += point.param_e * std::sin(point.x / (point.param_f * point.param_f + detail::eps)); 
  
            point.x *= weight;
            point.y *= weight;
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V15<data_t>::name {"waves"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V16 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V16(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.0, 1.0, 0.4})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            auto r_factor = (2.0 * weight) / (std::sqrt(point.x * point.x + point.y * point.y) + 1);
          
            //swap the x and y coordinates
            point.x = r_factor * point.y;
            point.y = r_factor * point.x;
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V16<data_t>::name {"fisheye"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V17 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V17(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.65, 0.0, 0.9})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            point.x += point.param_c * std::sin(std::tan(3 * point.y)); 
            point.y += point.param_f * std::sin(std::tan(3 * point.x)); 
  
            point.x *= weight;
            point.y *= weight;
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V17<data_t>::name {"popcorn"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V18 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V18(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.5, 1.0, 0.0})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            const data_t exp_factor = weight * std::exp(point.x - 1);
            const data_t dy = fflame_constants::PI * point.y;
            point.x = exp_factor * std::cos(dy);
            point.y = exp_factor * std::sin(dy);
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V18<data_t>::name {"exponential"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V19 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V19(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.0, 0.85, 0.5})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            const auto theta = std::atan2(point.x, point.y);
            const auto r_factor = weight * std::pow(std::sqrt(point.x * point.x + point.y * point.y), std::sin(theta));
          
            point.x = r_factor * std::cos(theta);
            point.y = r_factor * std::sin(theta);
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V19<data_t>::name {"power"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V20 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V20(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.75, 0.45, 0.2})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            point.x = weight * std::cos(fflame_constants::PI * point.x) * std::cosh(point.y);
            point.y = weight * -std::sin(fflame_constants::PI * point.x) * std::sinh(point.y);

            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V20<data_t>::name {"cosine"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V21 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V21(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.7, 0, 0.8})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            const auto theta = std::atan2(point.x, point.y);
            const auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
            const auto csq_param =  point.param_c*point.param_c + detail::eps;
            const auto coeff = weight * (std::fmod((r_factor + csq_param), (2*csq_param)) - csq_param + r_factor * (1 - csq_param));
  
            point.x = coeff * std::cos(theta);
            point.y = coeff * std::sin(theta);
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V21<data_t>::name {"rings"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V22 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V22(const double weight, const std::vector<data_t>& color = std::vector<data_t>{1.0, 0.2, 0.0})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            const auto r_factor = weight * std::sqrt(point.x * point.x + point.y * point.y);
            const auto theta = std::atan2(point.x, point.y);
          
            auto t_factor = fflame_constants::PI * (point.param_c*point.param_c + detail::eps);
            if(std::fmod((theta + point.param_f), t_factor) > (t_factor/2.0))
                t_factor *= -1;
  
            point.x = r_factor * (std::cos(theta + t_factor/2.0));
            point.y = r_factor * (std::sin(theta + t_factor/2.0));
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V22<data_t>::name {"fan"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V23 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V23(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0, 1.0, 0.5})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            const auto r_factor = (2.0 * weight) / (std::sqrt(point.x * point.x + point.y * point.y) + 1);
            point.x *= r_factor; 
            point.y *= r_factor; 
  
            detail::blend_colors(color, point);
        }
    };
    template <typename data_t>
    const std::string V23<data_t>::name {"eyefish"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V24 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V24(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.5, 0, 1.0})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            const auto r_factor = weight / ((point.x * point.x + point.y * point.y)/4.0 + 1);
            point.x *= r_factor; 
            point.y *= r_factor; 
  
            detail::blend_colors(color, point);        
        }
    };
    template <typename data_t>
    const std::string V24<data_t>::name {"bubble"};

    //////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct V25 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V25(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.75, 0.125, 0.45})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
            const auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
 
            point.x = r_factor * std::sin(point.x); 
            point.y *= r_factor; 
  
            detail::blend_colors(color, point);        
        }
    };
    template <typename data_t>
    const std::string V25<data_t>::name {"cylinder"};

    //////////////////////////////////////////////////////////////////
 
    template <typename data_t>
    struct V26 : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        V26 (const double weight, data_t blob_high, data_t blob_low, data_t blob_waves, const std::vector<data_t>& color = {0.2, 1.0, 0.1})
          : variant<data_t>(weight, color) 
        {
          blob_factor = (blob_high - blob_low)/2.0;
          blob_low = blob_low;
          blob_waves = blob_waves;
        }

        void apply_variant(flame_point<data_t>& point) override
        {
          const auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
          const auto theta = std::atan2(point.x, point.y);

          const auto coeff = weight * r_factor * (blob_low + blob_factor * (std::sin(blob_waves * theta) + 1));
          point.x = coeff * std::cos(theta);
          point.y = coeff * std::sin(theta);

          detail::blend_colors(color, point);
        }

        data_t blob_factor;
        data_t blob_low;
        data_t blob_waves;
    };
    template <typename data_t>
    const std::string V26<data_t>::name {"blob"};
   
/*
    //////////////////////////////////////////////////////////////////
    template <typename data_t>
    struct V : variant<data_t>
    {
        const static std::string name;
        using variant<data_t>::weight;
        using variant<data_t>::color;

        explicit V(const double weight, const std::vector<data_t>& color = std::vector<data_t>{0, 0, 0})
            : variant<data_t>(weight, color)
        {}

        void apply_variant(flame_point<data_t>& point) override
        {
        }
    };
    template <typename data_t>
    const std::string V<data_t>::name {};
*/

    ////////////////////////////////////////////////////////////
    //    this is the start of the parameteric equations
    ////////////////////////////////////////////////////////////


    //TODO still have tons more to go...

    template <template <class> class variant_base_t, template <class> class variant_t, 
          typename data_t, typename ... variant_args>
    variant_base_t<data_t>* generate_variant (variant_args&& ... args)
    {
        return new variant_t<data_t>(std::forward<variant_args>(args) ... );
    }

    template <typename data_t>
    struct variant_list
    {
        variant_list()
        {
            std::function<affine_fcns::variant<data_t>* ()> fobj = 
            []()
            {
                std::vector<data_t> default_color {0.25, 0.25, 1.0};
                return generate_variant<variant, V0, data_t>(1.0, default_color);
            };
            flame_maker.register_product(V0<data_t>::name, fobj);

            flame_maker.register_product(V1<data_t>::name, [](){return generate_variant<variant, V1, data_t>(1.0);});            
            flame_maker.register_product(V2<data_t>::name, [](){return generate_variant<variant, V2, data_t>(1.0);});            
            flame_maker.register_product(V3<data_t>::name, [](){return generate_variant<variant, V3, data_t>(1.0);});            
            flame_maker.register_product(V4<data_t>::name, [](){return generate_variant<variant, V4, data_t>(1.0);});            
            flame_maker.register_product(V5<data_t>::name, [](){return generate_variant<variant, V5, data_t>(1.0);});            
            flame_maker.register_product(V6<data_t>::name, [](){return generate_variant<variant, V6, data_t>(1.0);});            
            flame_maker.register_product(V7<data_t>::name, [](){return generate_variant<variant, V7, data_t>(1.0);});            
            flame_maker.register_product(V8<data_t>::name, [](){return generate_variant<variant, V8, data_t>(1.0);});            
            flame_maker.register_product(V9<data_t>::name, [](){return generate_variant<variant, V9, data_t>(1.0);});            

            flame_maker.register_product(V10<data_t>::name, [](){return generate_variant<variant, V10, data_t>(1.0);});            
            flame_maker.register_product(V11<data_t>::name, [](){return generate_variant<variant, V11, data_t>(1.0);});            
            flame_maker.register_product(V12<data_t>::name, [](){return generate_variant<variant, V12, data_t>(1.0);});            
            flame_maker.register_product(V13<data_t>::name, [](){return generate_variant<variant, V13, data_t>(1.0);});            
            flame_maker.register_product(V14<data_t>::name, [](){return generate_variant<variant, V14, data_t>(1.0);});            
            flame_maker.register_product(V15<data_t>::name, [](){return generate_variant<variant, V15, data_t>(1.0);});            
            flame_maker.register_product(V16<data_t>::name, [](){return generate_variant<variant, V16, data_t>(1.0);});            
            flame_maker.register_product(V17<data_t>::name, [](){return generate_variant<variant, V17, data_t>(1.0);});            
            flame_maker.register_product(V18<data_t>::name, [](){return generate_variant<variant, V18, data_t>(1.0);});            
            flame_maker.register_product(V19<data_t>::name, [](){return generate_variant<variant, V19, data_t>(1.0);});            

            flame_maker.register_product(V20<data_t>::name, [](){return generate_variant<variant, V20, data_t>(1.0);});            
            flame_maker.register_product(V21<data_t>::name, [](){return generate_variant<variant, V21, data_t>(1.0);});            
            flame_maker.register_product(V22<data_t>::name, [](){return generate_variant<variant, V22, data_t>(1.0);});            
            flame_maker.register_product(V23<data_t>::name, [](){return generate_variant<variant, V23, data_t>(1.0);});            
            flame_maker.register_product(V24<data_t>::name, [](){return generate_variant<variant, V24, data_t>(1.0);});            
            flame_maker.register_product(V25<data_t>::name, [](){return generate_variant<variant, V25, data_t>(1.0);});            
        }

        const static std::vector<std::string> variant_names;
        typedef Factory<affine_fcns::variant<data_t>, std::string, std::function<affine_fcns::variant<data_t>* ()>> flame_factory;
        flame_factory flame_maker;        
    };
    template <typename data_t>
    const std::vector<std::string> variant_list<data_t>::variant_names  
    {
        V0<data_t>::name,  V1<data_t>::name, V2<data_t>::name, V3<data_t>::name, V4<data_t>::name, 
        V5<data_t>::name,  V6<data_t>::name,  V7<data_t>::name,  V8<data_t>::name,  V9<data_t>::name, 
        V10<data_t>::name, V11<data_t>::name, V12<data_t>::name, V13<data_t>::name, V14<data_t>::name, 
        V15<data_t>::name, V16<data_t>::name, V17<data_t>::name, V18<data_t>::name, V19<data_t>::name,
        V20<data_t>::name, V21<data_t>::name, V22<data_t>::name, V23<data_t>::name, V24<data_t>::name,
        V25<data_t>::name 
    };

////////////////////////////////////////////////////////////////////////////////////////

    template <typename data_t>
    struct invoker
    {
        using flame_fcn = std::shared_ptr<variant<data_t>>; 
        invoker(std::vector<flame_fcn>&& fcn_default)
            : fcn(std::move(fcn_default))
        {
            flame_fcn_params<data_t> default_p {0, 0, 0, 0, 0, 0};

            for (int fcn_idx = 0; fcn_idx < fcn.size(); fcn_idx++) {
                affine_preparameters.push_back(default_p);
                affine_postparameters.push_back(default_p);    
            }
            /*
            for (auto f : fcn)
            {
                affine_preparameters.insert(std::make_pair(f, default_p));
                affine_postparameters.insert(std::make_pair(f, default_p));    
            }
            */
        }

        void invoke(const size_t fcn_idx, flame_point<data_t>& pt) const
        {
            //apply the selected function's parameters to the input point
            assert(fcn_idx < fcn.size());
            auto flame_function = fcn[fcn_idx];
            auto fflame_randeng = fflame_util::get_engine();

            //just for fun -- randomize the parameter lists too. Will want to change this to use the fast_rand 
            // -- apparently this takes ~15% of the TOTAL execution time (according to callgrind)
            static std::uniform_int_distribution<> fcn_dis(0, fcn.size()-1); 
            int pre_fcn_idx = fcn_dis(fflame_randeng);
            pt.apply_affine_params(affine_preparameters[pre_fcn_idx]);

            /*
            auto param_preit = affine_preparameters.find(fcn[fcn_dis(fflame_randeng)]);
            //auto param_preit = affine_preparameters.find(flame_function);
            if(param_preit != affine_preparameters.end())
                pt.apply_affine_params(param_preit->second);
            */

            flame_function->apply_variant (pt);

            //apply the post-processing affine transform parameters

            int post_fcn_idx = fcn_dis(fflame_randeng);
            pt.apply_affine_params(affine_postparameters[post_fcn_idx]);

            /*
            auto param_postit = affine_postparameters.find(fcn[fcn_dis(fflame_randeng)]);
            //auto param_postit = affine_postparameters.find(flame_function);
            if(param_postit != affine_postparameters.end())
                pt.apply_affine_params(param_postit->second);        
            */
        }

        //assign random values between [min_param, max_param) to all affine pre- and post-process parameters
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
                auto rng_engine = fflame_util::get_engine();
                flame_fcn_params<data_t> affine_preparams {
                    dis(rng_engine), dis(rng_engine), dis(rng_engine),
                    dis(rng_engine), dis(rng_engine), dis(rng_engine)
                };
                affine_preparameters [fcn_idx] = std::move(affine_preparams);

                /*
                auto param_preit = affine_preparameters.find(fcn[fcn_idx]);
                if(param_preit != affine_preparameters.end())
                {
                    flame_fcn_params<data_t> affine_params {
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine())
                    };
                    param_preit->second = std::move(affine_params);
                }
                */

                //make the post-parameters
                flame_fcn_params<data_t> affine_postparams {
                    dis(rng_engine), dis(rng_engine), dis(rng_engine),
                    dis(rng_engine), dis(rng_engine), dis(rng_engine)
                };
                affine_postparameters [fcn_idx] = std::move(affine_postparams);

                /*
                auto param_postit = affine_postparameters.find(fcn[fcn_idx]);
                if(param_postit != affine_postparameters.end())
                {
                    flame_fcn_params<data_t> affine_params {
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine()),
                        dis(fflame_util::get_engine())
                    };
                    param_postit->second = std::move(affine_params);
                }
                */                
            }
        }
        
        std::vector<flame_fcn> fcn;

        std::map<std::string, flame_fcn> fcn_finder;
        std::vector<data_t> fcn_probabilities;

        std::vector<flame_fcn_params<data_t>> affine_preparameters;
        std::vector<flame_fcn_params<data_t>> affine_postparameters;

        //std::map<flame_fcn, flame_fcn_params<data_t>> affine_preparameters;
        //std::map<flame_fcn, flame_fcn_params<data_t>> affine_postparameters;
    };
}

    /*
    //just a linear function
    template <typename data_t>
    void V0 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.25, 0.25, 1.0})
    {
        const static std::string name {"linear"};
        detail::blend_colors(color, point);

        point.x = weight * std::sin(point.x);
        point.y = weight * std::sin(point.y);
    }
    */


    /*
    template <typename data_t>
    void V1 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t>{1.0, 0.0, 0.0})
    {
        const static std::string name {"sinusoidal"};

        point.x = weight * std::sin(point.x);
        point.y = weight * std::sin(point.y);
        detail::blend_colors(color, point);
    }
    */

 /*
    template <typename data_t>
    void V2 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t>{0, 1.0, 0})
    {
        const static std::string name {"spherical"};
        
        auto r_factor = 1.0 / (point.x * point.x + point.y * point.y + detail::eps);
        point.x = weight * point.x * r_factor;
        point.y = weight * point.y * r_factor;
        detail::blend_colors(color, point);
    }
    */
   /*
    template <typename data_t>
    void V3 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.5, 0.1, 0.75})
    {
        const static std::string name {"swirl"};
        
        auto theta = point.x * point.x + point.y * point.y;
        point.x = weight * point.x * std::sin(theta) - point.y * std::cos(theta);
        point.y = weight * point.x * std::cos(theta) + point.y * std::sin(theta);
        detail::blend_colors(color, point);
    }
    */

    /*
    template <typename data_t>
    void V4 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t>{0, 0.75, 0.75})
    {
        const static std::string name {"horseshoe"};

        auto r_factor = 1.0 / (std::sqrt(point.x * point.x + point.y * point.y) + detail::eps);
     
        point.x = weight * (point.x - point.y)*(point.x + point.y) * r_factor;
        point.y = weight * 2 * point.x * point.y * r_factor; 
        detail::blend_colors(color, point);
    }
    */
    /*
    template <typename data_t>
    void V5 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t>{0.5, 0.5, 0.5})
    {   
        const static std::string name {"polar"};
    
        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);

        point.x = weight * theta / fflame_constants::PI;
        point.y = weight * r_factor - 1.0;

        detail::blend_colors(color, point);
    }
    */
    /*
    //handkerchief
    template <typename data_t>
    void V6 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.4, 0, 0.6})
    {   
        const static std::string name {"handkerchief"};
 
        auto r_factor = weight * std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);
       
        point.x = r_factor * std::sin(theta + r_factor);
        point.y = r_factor * std::cos(theta - r_factor);

        detail::blend_colors(color, point);
    }


    //heart 
    template <typename data_t>
    void V7 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0, 1.0, 0.3})
    {   
        const static std::string name {"heart"};
                
        auto r_factor = weight * std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);
       
        point.x = r_factor * std::sin(theta * r_factor);
        point.y = -r_factor * std::cos(theta * r_factor);       

        detail::blend_colors(color, point);
    }

    //disk
    template <typename data_t>
    void V8 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.9, 0.2, 0.2})
    {   
        const static std::string name {"disk"};
 
        auto factor = fflame_constants::PI * std::sqrt(point.x * point.x + point.y * point.y);
        auto coeff = weight * std::atan2(point.x * fflame_constants::PI, point.y * fflame_constants::PI) / fflame_constants::PI;

        point.x = coeff * std::sin(factor);
        point.y = coeff * std::cos(factor);       

        detail::blend_colors(color, point);
    }

    //spiral
    template <typename data_t>
    void V9 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.35, 0.25, 0.7})
    {   
        const static std::string name {"spiral"};
     
        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);

        point.x = weight * (std::cos(theta) + std::sin(r_factor)) / (r_factor + detail::eps);
        point.y = weight * (std::sin(theta) - std::cos(r_factor)) / (r_factor + detail::eps);
   
        detail::blend_colors(color, point);
    }

    //hyperbolic
    template <typename data_t>
    void V10 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.1, 1.0, 0.1})
    {   
        const static std::string name {"hyperbolic"};
 
        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);

        point.x = weight * std::sin(theta) / r_factor;
        point.y = weight * r_factor * std::cos(theta);
      
        detail::blend_colors(color, point);
    }

    //diamond
    template <typename data_t>
    void V11 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.15, 0.15, 1.0})
    {   
        const static std::string name {"diamond"};

        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);

        point.x = weight * std::sin(theta) * std::cos(r_factor);
        point.y = weight * std::cos(theta) * std::sin(r_factor);

        detail::blend_colors(color, point);
    }

    //Ex
    template <typename data_t>
    void V12 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.2, 0.4, 0.8})
    {   
        const static std::string name {"ex"};

        auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        auto theta = std::atan2(point.x, point.y);
        auto p0 = std::sin(theta + r_factor);
        auto p1 = std::cos(theta - r_factor);

        point.x = weight * r_factor * (p0*p0*p0 + p1*p1*p1);
        point.y = weight * r_factor * (p0*p0*p0 - p1*p1*p1);

        detail::blend_colors(color, point);
    }

    //julia
    template <typename data_t>
    void V13 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.8, 0.6, 0.2})
    {   
        const static std::string name {"julia"};

        //no, the double sqrt is not a typo
        auto r_factor = weight * std::pow(std::sqrt(std::sqrt(point.x * point.x + point.y * point.y)), 0.25);
        auto theta = std::atan2(point.x, point.y) / 2.0;

        //omega is "a random variable that's either 0 or pi"
        static std::uniform_int_distribution<> omega_dist(0, 1); 
        theta += fflame_constants::PI * omega_dist(fflame_util::get_engine());

        point.x = r_factor * std::cos(theta);
        point.y = r_factor * std::sin(theta);

        detail::blend_colors(color, point);
    }

    //bent
    template <typename data_t>
    void V14 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.25, 0.75, 0.25})
    {   
        const static std::string name {"bent"};
       
        if(point.y < 0)
            point.y /= 2;
        if(point.x < 0)
            point.x *= 2;

        point.x *= weight;
        point.y *= weight;

        detail::blend_colors(color, point);
    }

    //waves
    template <typename data_t>
    void V15 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.5, 0.625, 0.125})
    {   
        const static std::string name {"waves"};

        point.x += point.param_b * std::sin(point.y / (point.param_c * point.param_c + detail::eps));
        point.y += point.param_e * std::sin(point.x / (point.param_f * point.param_f + detail::eps)); 

        point.x *= weight;
        point.y *= weight;

        detail::blend_colors(color, point);
    }

    //fisheye
    template <typename data_t>
    void V16 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.1, 1.0, 0.4})
    {   
        const static std::string name {"fisheye"};

        auto r_factor = (2.0 * weight) / (std::sqrt(point.x * point.x + point.y * point.y) + 1);
        
        //swap the x and y coordinates
        point.x = r_factor * point.y;
        point.y = r_factor * point.x;

        detail::blend_colors(color, point);
    }

    //popcorn
    template <typename data_t>
    void V17 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.0, 0.65, 0.9})
    {   
        const static std::string name {"popcorn"};

        point.x += point.param_c * std::sin(std::tan(3 * point.y)); 
        point.y += point.param_f * std::sin(std::tan(3 * point.x)); 

        point.x *= weight;
        point.y *= weight;

        detail::blend_colors(color, point);
    }
   
    //exponential
    template <typename data_t>
    void V18 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.5, 0.9, 0.0})
    {   
        const static std::string name {"exponential"};

        const data_t exp_factor = weight * std::exp(point.x - 1);
        const data_t dy = fflame_constants::PI * point.y;
        point.x = exp_factor * std::cos(dy);
        point.y = exp_factor * std::sin(dy);

        detail::blend_colors(color, point);
    }
    
    //power
    template <typename data_t>
    void V19 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.0, 0.85, 0.45})
    {   
        const static std::string name {"power"};

        const auto theta = std::atan2(point.x, point.y);
        const auto r_factor = weight * std::pow(std::sqrt(point.x * point.x + point.y * point.y), std::sin(theta));
        
        point.x = r_factor * std::cos(theta);
        point.y = r_factor * std::sin(theta);

        detail::blend_colors(color, point);
    }
    
    //cosine
    template <typename data_t>
    void V20 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.45, 0.40, 0.25})
    {
        const static std::string name {"cosine"};

        point.x = weight * std::cos(fflame_constants::PI * point.x) * std::cosh(point.y);
        point.y = weight * -std::sin(fflame_constants::PI * point.x) * std::sinh(point.y);

        detail::blend_colors(color, point);
    }

    //rings
    template <typename data_t>
    void V21 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {0.7, 0, 0.8})
    {
        const static std::string name {"rings"};
        
        const auto theta = std::atan2(point.x, point.y);
        const auto r_factor = std::sqrt(point.x * point.x + point.y * point.y);
        const auto csq_param =  point.param_c*point.param_c + detail::eps;
        const auto coeff = weight * (((r_factor + csq_param) % (2*csq_param)) - csq_param + r_factor * (1 - csq_param))

        point.x = coeff * std::cos(theta);
        point.y = coeff * std::sin(theta);

        detail::blend_colors(color, point);
    }

    //fan
    template <typename data_t>
    void V22 (flame_point<data_t>& point, const double weight, const std::vector<data_t>& color = std::vector<data_t> {1.0, 0.2, 0})
    {
        const static std::string name {"fan"};

        const auto r_factor = weight * std::sqrt(point.x * point.x + point.y * point.y);
        const auto theta = std::atan2(point.x, point.y);
        
        auto t_factor = fflame_constants::PI * (point.param_c*point.param_c + detail::eps);
        if((theta + point.param_f) % t_factor > (t_factor/2.0))
            t_factor *= -1;

        point.x = r_factor * (std::cos(theta + t_factor/2.0));
        point.y = r_factor * (std::sin(theta + t_factor/2.0));

        detail::blend_colors(color, point);
    }
      */
    /*
    template <typename fcn_t, typename ... fcn_params>
    std::function<void(flame_point<data_t>& point)> make_variant(fcn_t, fcn_params&& ... args)
    {
        std::function<void(flame_point<data_t>& point)> fcn = [&](flame_point<data_t>& point)
        {
            return fcn_t(std::ref(point), std::forward<fcn_params>(args) ...);
        }
    }*/

    ////////////////////////////////////////////////////////////////////////////////////////
    //will have to keep this updated as we add new variants
    //constexpr uint8_t num_variants = 6; 


#endif
