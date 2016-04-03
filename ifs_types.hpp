/* ifs_types.hpp -- part of the fractal flames implementation 
 *
 * Copyright (C) 2015 Alrik Firl 
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef FF_IFS_CONSTANTS_HPP
#define FF_IFS_CONSTANTS_HPP

#include "util/ff_utils.hpp"

#include <random>
#include <tuple>
#include <thread>
#include <memory>
#include <iostream>
#include <map>
#include <type_traits>

namespace fflame_constants
{
    constexpr int num_samples = 1000;
    constexpr int max_iter = 20000;

    //generate real points within [-1, 1) to use for seeding
    constexpr double min_pt = -1.0;
    constexpr double max_pt = 1.0;

    constexpr int imheight = 1024;
    constexpr int imwidth = 1024;
    constexpr double gamma = 2.2; //1.2;
    constexpr double gamma_factor = 1.0/gamma;

    constexpr double PI = 3.14159265358979323; 
}

/*
template <class VariantType,
           typename KeyType,
           typename VariantCreator> 
class flame_factory
{
public:
    VariantType* create_variant(const KeyType& id)
    {
        auto variant_iter = creator_map.find(id);
        if(variant_iter != creator_map.end()) {
            return (variant_iter->second)();
        }
        return nullptr;
    }

    bool register_variant(const KeyType& id, VariantCreator creator)
    {
        return creator_map.insert(std::pair<KeyType, VariantCreator>(id, creator)).second;
    }

    bool unregister_variant(const KeyType& id)
    {
        return (creator_map.erase(id) == 1);
    }

private:
    std::map<KeyType, VariantCreator> creator_map;
};
*/

//assume that we won't need to manage any of the resources here (i.e. all data is stack-allocated
//and/or managed by other entities)
template <typename pixel_t>
class flame_frame
{
public:
    using T = pixel_t;

    flame_frame()
        : rows(0), cols(0), data(nullptr), manage_data(true)
    {}

    flame_frame(int height, int width)
        : rows(height), cols(width), manage_data(true)
    {
        data = new pixel_t [rows * cols];
    }
    flame_frame(int height, int width, const pixel_t& initializer)
        : rows(height), cols(width), manage_data(true)
    {
        data = new pixel_t [rows * cols];
        if(std::is_pod<pixel_t>::value) {
            std::fill(data, data + rows * cols, initializer);
        } else {
            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    data [r * cols + c] = initializer;
                }
            }
        }
    }

    flame_frame(int height, int width, pixel_t* data_p)
        : rows(height), cols(width), data(data_p), manage_data(false)
    {}


    flame_frame(const flame_frame& other) 
        : rows(other.rows), cols(other.cols), data(other.data), manage_data(other.manage_data)
    {}

    flame_frame(flame_frame&& other)
        : flame_frame() 
    {
        swap(*this, other);
    }

    flame_frame& operator=(flame_frame other)
    {
        swap(*this, other);
        return *this;
    }

    ~flame_frame()
    {
        if(manage_data) {
            delete [] data;
        }
    }

    inline pixel_t* ptr (const int row)
    {
        return data + row * cols;
    }

    inline const pixel_t* ptr (const int row) const
    {
        return data + row * cols;
    }

    inline pixel_t& at (const int row, const int col)
    {
        return *(data + row * cols + col);

    }

    inline const pixel_t& at (const int row, const int col) const
    {
        return *(data + row * cols + col);
    }

    friend void swap(flame_frame& lhs, flame_frame& rhs)
    {
        using std::swap; 
        swap(lhs.rows, rhs.rows); 
        swap(lhs.cols, rhs.cols);
        swap(lhs.data, rhs.data);
        swap(lhs.manage_data, rhs.manage_data);
    }

    int rows;
    int cols;
    pixel_t* data;

private:
    bool manage_data;
};

template <typename pixel_t>
struct histogram_info
{
    histogram_info()
        : color{0, 0, 0}
    {
        frequency_count = 0;
    }

    histogram_info(pixel_t px_info, int freq_val)
        : color(px_info)
    {
        frequency_count = freq_val;
    }

    void update(const histogram_info& other)
    {
        color += other.color; 
        frequency_count += other.frequency_count;
    }

    void update(const pixel_t& px_info, const int freq_info)
    {
        color += px_info; 
        frequency_count += freq_info;
    }

    void reset()
    {
        color[0] = 0;
        color[1] = 0;
        color[2] = 0;
        frequency_count = 0;
    }

    pixel_t color;
    int frequency_count;
};

template <typename data_t>
struct flame_fcn_params
{
    flame_fcn_params(data_t a = 0, data_t b = 0, data_t c = 0, data_t d = 0, data_t e = 0, data_t f = 0)
        : x_affine(std::make_tuple(a,b,c)), y_affine(std::make_tuple(d,e,f))
    {}

    std::tuple<data_t, data_t, data_t> x_affine;
    std::tuple<data_t, data_t, data_t> y_affine;
};

//eventually will want to try this as spherical coordinates?
template <typename data_t>
struct flame_point
{
    flame_point(const data_t row, const data_t col)
      : y(row), x(col), color {0.0, 0.0, 0.0},
        param_a(0.0), param_b(0.0), param_c(0.0), 
        param_d(0.0), param_e(0.0), param_f(0.0)
    {}

    void apply_affine_params(const flame_fcn_params<data_t>& affine_params)
    {
        param_a = std::get<0>(affine_params.x_affine); 
        param_b = std::get<1>(affine_params.x_affine); 
        param_c = std::get<2>(affine_params.x_affine); 
        param_d = std::get<0>(affine_params.y_affine); 
        param_e = std::get<1>(affine_params.y_affine);
        param_f = std::get<2>(affine_params.y_affine);

        x = param_a * x + param_b * y + param_c; 
        y = param_d * x + param_e * y + param_f;     
    }

    data_t y;
    data_t x;
    //arranged as r, g, b
    data_t color [3];

    //also store the last-applied affine variation parameters for the dependant variations
    data_t param_a, param_b, param_c, param_d, param_e, param_f;
};

//helper threading class for running the fractal flames.
//carries the necessary per-thread state (since I don't have
//thread_local supported on my compiler apparently)
class flame_thread
{
public:
    flame_thread(uint64_t seed0, uint64_t seed1)
       : rand_gen(seed0, seed1), fthread(nullptr)
    {}

    flame_thread(const flame_thread&) = delete;
    flame_thread& operator= (const flame_thread&) = delete;

    
    flame_thread(flame_thread&& other)
        : rand_gen(other.rand_gen)
    {
        fthread = std::move(other.fthread);

        //what to do if the current object's thread object is running?
        other.rand_gen = fflame_util::fast_rand(0, 0);
        other.fthread = nullptr;
    }
    
    ~flame_thread()
    {
        //note: it's not a very good idea to do it this way. We would have to interrupt the thread
        //first. Plus I'm pretty sure join throws exceptions. It'll be fine for the current usage,
        //but not in any other
        if(fthread)
        {
            std::cout << "DTOR Finishing thread " << fthread->get_id() << std::endl;
            fthread->join();
        }
    }

    template <typename fcn_t, typename ... fcnargs_t>
    std::tuple<bool, std::thread::id> do_flame(fcn_t fcn, fcnargs_t&& ... args)
    {
        //check if the thread already exists
        if(fthread) {
            return std::make_tuple(false, fthread->get_id());
        }

        fthread = std::unique_ptr<std::thread>(new std::thread(fcn, std::forward<fcnargs_t>(args)..., std::ref(rand_gen)));
        return std::make_tuple((fthread != nullptr), fthread->get_id());
    }

    void finish_flame()
    {
        std::cout << "Finishing thread " << fthread->get_id() << std::endl;
        fthread->join();    
        fthread.reset(nullptr);
    }

private:    
    fflame_util::fast_rand rand_gen;
    std::unique_ptr<std::thread> fthread;
};

#endif
