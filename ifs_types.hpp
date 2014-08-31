#ifndef FF_IFS_CONSTANTS_HPP
#define FF_IFS_CONSTANTS_HPP

#include <random>
#include <tuple>
#include <thread>
#include <memory>
#include <iostream>

namespace fflame_constants
{
    constexpr int num_samples = 1000;
    constexpr int max_iter = 20000;

    //generate real points within [-1, 1) to use for seeding
    constexpr double min_pt = -1.0;
    constexpr double max_pt = 1.0;

    constexpr int imheight = 1024;
    constexpr int imwidth = 1024;

    //not sure if having differently-sized histograms is useful, but for now...
    constexpr int hist_height = imheight;
    constexpr int hist_width = imwidth;

    constexpr double gamma = 2.2;
    constexpr double gamma_factor = 1.0/gamma;

    constexpr double PI = 3.14159265358979323; 
}

namespace fflame_randutil
{
    std::default_random_engine& get_engine()
    {
        static std::random_device rdev{};
        static std::default_random_engine eng{rdev()};
        return eng;
    }

    //we assume the min value in the range is 0
    struct fast_rand
    {
        fast_rand(uint64_t seed0, uint64_t seed1)
            : s{seed0, seed1}
        {}

        uint64_t xorshift128plus(uint64_t max_val) {
            uint64_t s1 = s[ 0 ];
            const uint64_t s0 = s[ 1 ];
            s[ 0 ] = s0;
            s1 ^= s1 << 23;
            return (( s[ 1 ] = ( s1 ^ s0 ^ ( s1 >> 17 ) ^ ( s0 >> 26 ) ) ) + s0) % max_val;
        }
        uint64_t s[2];
    };

    /*
    static int fcn_rand(const int min, const int max) {
        static std::thread_local std::mt19937 generator;
        std::uniform_int_distribution<int> distribution(min,max);
        return distribution(generator);
    }
    */

}

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

    void update(const pixel_t& px_info, const int freq_info)
    {
        color += px_info; 
        frequency_count += freq_info;
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
private:    
    fflame_randutil::fast_rand rand_gen;
    std::unique_ptr<std::thread> fthread;

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
        other.rand_gen = fflame_randutil::fast_rand(0, 0);
        other.fthread = nullptr;
    }
    

    ~flame_thread()
    {
        //note: it's not a very good idea to do it this way. We would have to interrupt the thread
        //first. Plus I'm pretty sure join throws exceptions. It'll be fine for the current usage,
        //but not in any other
        if(fthread)
        {
            std::cout << "Finishing thread " << fthread->get_id() << std::endl;
            fthread->join();
        }
    }

    template <typename fcn_t, typename ... fcnargs_t>
    bool do_flame(fcn_t fcn, fcnargs_t&& ... args)
    {
        //check if the thread already exists
        if(fthread)
            return false;

        fthread = std::unique_ptr<std::thread>(new std::thread(fcn, std::forward<fcnargs_t>(args)..., std::ref(rand_gen)));
        return (fthread != nullptr);
    }

    void finish_flame()
    {
        std::cout << "Finishing thread " << fthread->get_id() << std::endl;
        fthread->join();    
        fthread.reset(nullptr);
    }
};

#endif
