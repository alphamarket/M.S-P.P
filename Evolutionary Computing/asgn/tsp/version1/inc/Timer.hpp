#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <type_traits>

namespace TSPGA {
    template<class T = std::chrono::milliseconds>
    class Timer
    {
    public:
        Timer() : beg_(clock_::now()) {
            static_assert(std::is_class<std::chrono::duration<int64_t>>::value, "The <T> should be a type of `std::chrono::duration<int64_t>`");
        }
        void reset() { beg_ = clock_::now(); }
        double elapsed() const {
            return std::chrono::duration_cast<T>
                (clock_::now() - beg_).count(); }

    private:
        typedef std::chrono::high_resolution_clock clock_;
        std::chrono::time_point<clock_> beg_;
    };
    typedef Timer<std::chrono::nanoseconds>     Timer_nan;
    typedef Timer<std::chrono::microseconds>    Timer_mic;
    typedef Timer<std::chrono::milliseconds>    Timer_mil;
    typedef Timer<std::chrono::seconds>         Timer_sec;
    typedef Timer<std::chrono::minutes>         Timer_min;
    typedef Timer<std::chrono::hours>           Timer_hrs;
}

#endif // TIMER_H
