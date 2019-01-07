#ifndef _SOFTLIB_PROGRESS_TRACKER_H
#define _SOFTLIB_PROGRESS_TRACKER_H

#include <chrono>
#include <string>
#include <softlib/config.h>
#include <softlib/Timer.h>

class ProgressTracker {
    public:
        enum ProgressType {
            PROGRESS_LINES
        };
        enum TimeUnit {
            DAYS=1,
            HOURS=2,
            MINUTES=3,
            SECONDS=4,
            MILLISECONDS=5
        };
    private:
        Timer *start, *last;

        bool colors = false;    // If true, output is formatted with fancy colors
        bool estimate_remaining = false;    // Estimate how much time remains
        long long int
            dprint,          // Number of counts per progress point
            nprint,          // Number of progress points to print
            modlevel,        // If (total%nprint > 0), we add 1 when
                             // incrementing 'next' until this value
                             // has been reached
            total,           // Total number of progress points
            next,            // Next point at which to emit progress
            nextindex;       // Index (from 1 to nprint) of next point

        // Running times of three previous blocks
        //   dt1 = t_{i-1} - t_{i-2}
        //   dt2 = t_{i-2} - t_{i-3}
        long long int dt1, dt2;

        unsigned int ndigits;   // Number of digits in 'nprint'

        ProgressType type = PROGRESS_LINES;
    public:
        ProgressTracker(
            const long long int, const long long int, const ProgressType,
            const bool colors=false, const bool estimate=false
        );
        ~ProgressTracker();

        void PrintProgress(const long long int);
        void PrintProgressLines(const long long int);

        void UpdateNext();

        long long int DiffPassedTime();
        long long int EstimateRemainingTimeAverage(const long long int);
        long long int EstimateRemainingTimeDeterministic(const long long int);
        long long int TotalPassedTime();
        long long int PassedTime(Timer*);

        std::string TimeString(const long long int, const bool color=false, const TimeUnit tu=TimeUnit::HOURS);
        std::string FormattedTimeString(const long long int, const long long int, const long long int, const long long int, const long long int, const bool color=false);
};

#endif/*_SOFTLIB_PROGRESS_TRACKER_H*/
