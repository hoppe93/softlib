#ifndef _TIMER_H
#define _TIMER_H

#include <chrono>
#include <string>
#include <softlib/config.h>


class Timer {
    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> tic;
		std::chrono::time_point<std::chrono::high_resolution_clock> toc;

		std::chrono::high_resolution_clock::duration offset=std::chrono::high_resolution_clock::duration::zero();

		bool stopped=false;
    public:
        Timer();

        struct time {
            std::chrono::hours days;
            std::chrono::hours hours;
            std::chrono::minutes minutes;
            std::chrono::seconds seconds;
            std::chrono::milliseconds milliseconds;
            std::chrono::microseconds microseconds;
            std::chrono::seconds total_sec;
            std::chrono::microseconds tp;
        };

		void Continue();
        slibreal_t GetMilliseconds() const;
        static slibreal_t GetMilliseconds(struct time&);
        slibreal_t GetMicroseconds() const;
        static slibreal_t GetMicroseconds(struct time&);
        struct time GetTimeStruct() const;
        void Reset();
		void Stop();

        std::string ToString() const;
        static std::string ToString(struct time&);

        std::chrono::time_point<std::chrono::high_resolution_clock>& GetClock() { return tic; }
};

#endif/*_TIMER_H*/
