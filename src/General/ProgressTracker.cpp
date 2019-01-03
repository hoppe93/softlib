/**
 * :: A PROGRESS TRACKER/PRINTER
 *
 * The particle generator can print progress in two different ways:
 *   
 *   (i)  Using a progress bar (on terminals that support '\b' and '\r')
 *   (ii) Using continuous output to 'stdout' (aka SOFTv1 style).
 *
 * Mode (i) is automatically enabled when compiling with the macro
 * 'COLOR_TERMINAL' defined.
 */

#include <cstdio>
#include <chrono>
#include <sstream>
#include <string>
#include <softlib/config.h>
#include <softlib/ProgressTracker.h>

using namespace std;

/**
 * Initialize the progress tracker.
 *
 * n:      Total number of steps.
 * nprint: Number of times to print progress.
 * type:   Type of progress to output.
 * colors: If true, enables colored terminal output.
 */
ProgressTracker::ProgressTracker(
    const size_t n, const size_t nprint, const ProgressType type,
    const bool colors, const bool estimate
) {
    this->total = n;
    this->nprint = nprint;
    this->type = type;
    this->colors = colors;
    this->estimate_remaining = estimate;

    this->dprint = n / nprint;
    this->modlevel = (n % nprint) * (this->dprint+1);

    this->start = new Timer();
    this->last = new Timer();

    this->next = 0;
    this->nextindex = 0;
    UpdateNext();
}

/**
 * Destructor.
 */
ProgressTracker::~ProgressTracker() {
    delete this->last;
    delete this->start;
}

/**
 * Update the value of the next point to print progress at.
 */
void ProgressTracker::UpdateNext() {
    if (this->next < this->modlevel)
        this->next += this->dprint + 1;
    else
        this->next += this->dprint;

    this->nextindex++;
}

/**
 * Print progress (if necessary).
 */
void ProgressTracker::PrintProgress(const size_t i) {
    if (i >= this->next) {
        if (this->type == PROGRESS_LINES)
            PrintProgressLines();

        this->UpdateNext();
    }
}

/**
 */
void ProgressTracker::PrintProgressLines() {
    string remaining;
    string timestring = TimeString(TotalPassedTime());

    if (estimate_remaining)
        remaining = TimeString(EstimateRemainingTime());

    if (colors) {
        if (estimate_remaining)
            printf(
                "\x1B[1;32mPROGRESS\x1B[0m (%zu/%zu) :: %s passed, estimated %s remaining\n",
                this->nextindex, this->nprint,
                timestring.c_str(), remaining.c_str()
            );
        else
            printf(
                "\x1B[1;32mPROGRESS\x1B[0m (%zu/%zu) :: %s passed\n",
                this->nextindex, this->nprint,
                timestring.c_str()
            );
    } else {
        if (estimate_remaining)
            printf(
                "PROGRESS (%zu/%zu) :: %s passed, estimated %s remaining\n",
                this->nextindex, this->nprint,
                timestring.c_str(), remaining.c_str()
            );
        else
            printf(
                "PROGRESS (%zu/%zu) :: %s passed\n",
                this->nextindex, this->nprint,
                timestring.c_str()
            );
    }
}

/**
 * Estimate time remaining (assuming equal time
 * per progress point).
 */
size_t ProgressTracker::EstimateRemainingTime() {
    size_t passed = TotalPassedTime();
    size_t rem = (size_t)(((slibreal_t)passed) * (((slibreal_t)this->nprint)/((slibreal_t)this->nextindex)));

    return rem;
}

/**
 * Get the total time passed since this
 * object was initialized.
 */
size_t ProgressTracker::TotalPassedTime() { return PassedTime(this->start); }
/**
 * Get the amount of time that passed since
 * the last progress point was reported.
 */
size_t ProgressTracker::DiffPassedTime()  { return PassedTime(this->last); }

/**
 * Compute the amount of time that has passed
 * since the given Timer object was last reset.
 *
 * time: Timer object representing the start time.
 */
size_t ProgressTracker::PassedTime(Timer *tim) {
    chrono::time_point<chrono::high_resolution_clock> clk = chrono::system_clock::now();
    chrono::seconds sec = chrono::duration_cast<chrono::seconds>(clk - tim->GetClock());

    return ((size_t)sec.count());
}

/**
 * Convert a given number of seconds to a string with
 * the format
 *
 *   DdHHhMMmSSs
 *
 * Where
 *
 *   D  = Number of days
 *   HH = Number of hours (with preceeding 0)
 *   MM = Number of minutes (with preceeding 0)
 *   SS = Number of seconds (with preceeding 0)
 *
 * Zero times are omitted, unless a time with greater period
 * is non-zero (i.e. if D is returned, HH, MM and SS will always
 * be returned as well). If all values are zero, then "0s" is
 * returned.
 */
string ProgressTracker::TimeString(const size_t sec) {
    /*chrono::hours days   = chrono::duration_cast<hours>(sec) / 24;
    chrono::hours hours  = chrono::duration_cast<hours>(sec - days*24);
    chrono::minutes mins = chrono::duration_cast<minutes>(sec - days*24 - hours);
    chrono::seconds secs = chrono::duration_cast<seconds>(sec - days*24 - hours - mins);*/

    const slibreal_t
        MINUTE = 60.0,
        HOUR   = 60.0*MINUTE;
        //DAY    = 24.0*HOUR;

    //size_t days  = (size_t)(sec/DAY);
    //size_t hours = (size_t)((sec - days*DAY)/HOUR);
    //size_t mins  = (size_t)((sec - days*DAY - hours*HOUR)/MINUTE);
    //size_t secs  = (size_t) (sec - days*DAY - hours*HOUR - mins*MINUTE);
    size_t hours = (size_t)(sec/HOUR);
    size_t mins  = (size_t)((sec - hours*HOUR)/MINUTE);
    size_t secs  = (size_t) (sec - hours*HOUR - mins*MINUTE);

    ostringstream os;
    bool force = false;

    /*if (days > 0) {
        os << days << "d";
        force = true;
    }*/

    if (force || hours > 0) {
        if (hours < 10 && force)
            os << "0";

        os << hours << "h";
        force = true;
    }
    
    // MINUTES
    if (force || mins > 0) {
        if (mins < 10 && force)
            os << "0";

        os << mins << "m";
    }

    // SECONDS
    if (secs < 10 && force)
        os << "0";
    
    os << secs << "s";

    return os.str();
}

