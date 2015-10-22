/* ff_utils.hpp -- part of the fractal flames implementation 
 *
 * Copyright (C) 2015 Alrik Firl 
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef FRACTAL_FLAMES_UTILS_HPP
#define FRACTAL_FLAMES_UTILS_HPP

#include <random>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <iostream>

namespace fflame_util
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
        s[0] = s0;
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

class semaphore
{
public:    
    semaphore(int count = 0)
        : sem_count(count)
    {}

    void signal()
    {
        std::unique_lock<std::mutex> s_lock (sem_mtx);
        sem_count++;
        sem_cv.notify_one();
    }

    void wait()
    {
        std::unique_lock<std::mutex> s_lock (sem_mtx);
        sem_cv.wait(s_lock, [this]() 
                {
                    return sem_count > 0;
                });
        sem_count--;
    }
private:    
    std::mutex sem_mtx;
    std::condition_variable sem_cv;
    int sem_count;
};


class Barrier
{
private:
    std::mutex _mutex;
    std::condition_variable _cv;
    std::size_t _count;
public:
    explicit Barrier(std::size_t count) : _count{count} { }
    void Wait()
    {
        std::unique_lock<std::mutex> lock{_mutex};
        if (--_count == 0) {
            _cv.notify_all();
        } else {
            _cv.wait(lock, [this] { return _count == 0; });
        }
    }
};

class barrier
{
public:
    barrier(int count)
        : barrier_count(count)
    {}

    void wait()
    {
        std::unique_lock<std::mutex> b_lock (barrier_mtx);

        barrier_count--;
		if(barrier_count <= 0) {
            barrier_cv.notify_all();
		} else {
            barrier_cv.wait(b_lock, [this]()
            {
                return barrier_count <= 0;
            });
		}
    }

    void reset(int count)
    {
        std::unique_lock<std::mutex> b_lock (barrier_mtx);
        barrier_count = count;
    }

private:
    std::mutex barrier_mtx;
    std::condition_variable barrier_cv;
    int barrier_count;
};

} //namespace fflame_randutil
#endif
