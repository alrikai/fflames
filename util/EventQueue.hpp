#pragma once

#include <stdio.h>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include <atomic>
#include <thread>

template <typename EventType>
class EventQueue
{
public:
    explicit EventQueue(const int max_sz = 500, const int max_wait = 0) 
        : buffer_size(max_sz), timeout_len(max_wait) 
    {}

    ~EventQueue()
    {
        //clean out the frames in the Queue (NOTE: do we need this, since the elements are unique_ptrs?)
        std::lock_guard<std::mutex> lock(cleanup_);
        bool has_data = empty();
        while(has_data)
        {
            //since the queue returns a unique_ptr, this will invoke the element's destructor automatically
            pop(has_data);
        }
    }

    EventQueue(const EventQueue&) = delete;
    EventQueue& operator= (const EventQueue<EventType>&) = delete;
    
    void push(std::unique_ptr<EventType> data);
    std::unique_ptr<EventType> pop(bool& got_data);
    bool empty(void);
    
    bool has_started()
    {
        return queue_started.load(std::memory_order_seq_cst);
    }

    EventQueue(EventQueue&& other) 
        : EventQueue(nullptr)
    {
        std::lock_guard<std::mutex> lock(other.dlock_);
        
        swap(*this, other);
    }

    friend void swap(EventQueue& first, EventQueue& second)
    {
        using std::swap;

        swap(first.buffer_, second.buffer_);
        swap(first.buffer_size, second.buffer_size);
        swap(first.timeout_len, second.timeout_len); 
        swap(first.queue_started, second.queue_started); 
    }

private:
    std::queue<std::unique_ptr<EventType>> buffer_;
    mutable std::mutex dlock_;
    std::condition_variable dcond_;

    mutable std::mutex cleanup_;
    //maximum buffer size, will circle back around if maximum size reached
    int buffer_size;
    int timeout_len;

    std::atomic<int> num_producers;
    std::atomic<int> num_consumers;
    std::atomic<bool> queue_started;
};

template <typename EventType>
void EventQueue<EventType>::push(std::unique_ptr<EventType> data)
{
    std::lock_guard<std::mutex> lock(dlock_);
    //delete the oldest element to be replaced 
    if(buffer_.size() >= buffer_size)
    {
        std::cout << "Queue Full, Deleting Oldest Frame -- Thread " << std::this_thread::get_id() << std::endl;
        buffer_.pop();
    }

    buffer_.push(std::move(data));
    dcond_.notify_one();
}

template <typename EventType>
std::unique_ptr<EventType> EventQueue<EventType>::pop(bool& got_data)
{
    std::unique_lock<std::mutex> lock(dlock_);

    //to avoid unnecessary timeouts (even if it's short)
    if(buffer_.empty())
    {
        got_data = false;
        return nullptr;
    }

    //have a bounded wait, such that the thread waits for (up to) the timeout_len for data
    if(dcond_.wait_for(lock, std::chrono::milliseconds(timeout_len), [this] 
                { return !buffer_.empty(); })) 
    {
        //move the data element from the buffer 
        std::unique_ptr<EventType> data = std::move(buffer_.front());
        //delete the (now empty) object at the front of the queue
        buffer_.pop();
        got_data = true;
        return data;
    }
    else
    {
        got_data = false;
        return nullptr;
    }
}

template <typename EventType>
bool EventQueue<EventType>::empty()
{ 
    std::lock_guard<std::mutex> lock(dlock_);
    return buffer_.empty(); 
}
