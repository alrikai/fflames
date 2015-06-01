#ifndef FRACTAL_FLAMES_UTILS_HPP
#define FRACTAL_FLAMES_UTILS_HPP

#include <random>
#include <thread>
#include <mutex>
#include <condition_variable>

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
        barrier_cv.wait(b_lock, [this]()
            {
                return barrier_count <= 0;
            });
        barrier_cv.notify_all();
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
