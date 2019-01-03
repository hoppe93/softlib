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

    private:
        Timer *start, *last;

        bool colors = false;    // If true, output is formatted with fancy colors
        bool estimate_remaining = false;    // Estimate how much time remains
        size_t dprint,          // Number of counts per progress point
               nprint,          // Number of progress points to print
               modlevel,        // If (total%nprint > 0), we add 1 when
                                // incrementing 'next' until this value
                                // has been reached
               total,           // Total number of progress points
               next,            // Next point at which to emit progress
               nextindex;       // Index (from 1 to nprint) of next point

        unsigned int ndigits;   // Number of digits in 'nprint'

        ProgressType type = PROGRESS_LINES;
    public:
        ProgressTracker(
            const size_t, const size_t, const ProgressType,
            const bool colors=false, const bool estimate=false
        );
        ~ProgressTracker();

        void PrintProgress(const size_t);

        void PrintProgressBar();
        void PrintProgressLines();

        void UpdateNext();

        size_t DiffPassedTime();
        size_t EstimateRemainingTime();
        size_t TotalPassedTime();
        size_t PassedTime(Timer*);

        std::string TimeString(const size_t);
};

#endif/*_SOFTLIB_PROGRESS_TRACKER_H*/
