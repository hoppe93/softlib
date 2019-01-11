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

#include <chrono>
#include <cmath>
#include <cstdio>
#include <omp.h>
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
    const long long int n, const long long int nprint,
    const ProgressType type, const bool colors, const bool estimate
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
void ProgressTracker::PrintProgress(const long long int i) {
    const long long int index = this->nextindex;
    bool returnImmediately = false;

    #pragma omp critical (ProgressTracker_Print)
    {
        if (i < this->next)
            returnImmediately = true;
        else
            this->UpdateNext();
    }

    if (returnImmediately)
        return;

    if (this->type == PROGRESS_LINES)
        PrintProgressLines(index);

    this->last->Reset();
}

/**
 */
void ProgressTracker::PrintProgressLines(const long long int index) {
    string remaining, totalrem;
    long long int totalpassed = TotalPassedTime();
    string absolute = TimeString(totalpassed, colors);
    string relative = TimeString(DiffPassedTime(), colors);

    if (estimate_remaining) {
        long long int est = EstimateRemainingTimeDeterministic(index);
        //long long int est = EstimateRemainingTimeAverage(index);
        remaining = TimeString(est, colors);
        totalrem  = TimeString(est+totalpassed, colors);
    }

    if (colors) {
        if (estimate_remaining)
            printf(
                "\x1B[1;32mPROGRESS\x1B[0m (%lld/%lld) :: %s / %s passed, %s remaining (%s total)\n",
                index, this->nprint,
                absolute.c_str(), relative.c_str(), remaining.c_str(), totalrem.c_str()
            );
        else
            printf(
                "\x1B[1;32mPROGRESS\x1B[0m (%lld/%lld) :: %s / %s passed\n",
                index, this->nprint,
                absolute.c_str(), relative.c_str()
            );
    } else {
        if (estimate_remaining)
            printf(
                "PROGRESS (%lld/%lld) :: %s / %s passed, %s remaining (%s total)\n",
                index, this->nprint,
                absolute.c_str(), relative.c_str(), remaining.c_str(), totalrem.c_str()
            );
        else
            printf(
                "PROGRESS (%lld/%lld) :: %s / %s passed\n",
                index, this->nprint,
                absolute.c_str(), relative.c_str()
            );
    }
}

/**
 * Estimate time remaining after executing operation by
 * averaging over previous durations.
 */
long long int ProgressTracker::EstimateRemainingTimeAverage(const long long int i) {
    long long int ti = TotalPassedTime();
    long long int N = this->nprint;

    return (((N-i)*ti)/i);
}

/**
 * Estimate time remaining after executing operation i,
 * assuming the duration of each operation varies deterministically.
 * Returns the number of remaining milliseconds.
 */
long long int ProgressTracker::EstimateRemainingTimeDeterministic(const long long int i) {
    long long int dt = DiffPassedTime(), d2t=0, d3t=0, eT;
    long long int Ni = this->nprint - i;
    long long int NNp12, NNp1Np26;

    eT = Ni * dt;

    if (i >= 2) {
        d2t = dt - this->dt1;

        if (Ni%2==0) NNp12 = (Ni/2)*(Ni+1);
        else NNp12 = Ni*((Ni+1)/2);

        if (d2t*dt > 0)
            eT += NNp12 * d2t;
    }

    /*if (i >= 3) {
        d3t = (dt + this->dt2 - 2*this->dt1);

        if (Ni%3==0)
            NNp1Np26 = (Ni/3)*((Ni+1)/2)*(Ni+2);
        else if ((Ni+1)%3==0)
            NNp1Np26 = (Ni/2)*((Ni+1)/3)*(Ni+2);
        else
            NNp1Np26 = Ni*((Ni+1)/2)*((Ni+2)/3);

        if (d2t*dt > 0)
            eT += NNp1Np26 * d3t;
    }*/

    this->dt2 = this->dt1;
    this->dt1 = dt;

    return eT;
}

/**
 * Get the total time passed since this
 * object was initialized.
 */
long long int ProgressTracker::TotalPassedTime() { return PassedTime(this->start); }
/**
 * Get the amount of time that passed since
 * the last progress point was reported.
 */
long long int ProgressTracker::DiffPassedTime()  { return PassedTime(this->last); }

/**
 * Compute the amount of time that has passed
 * since the given Timer object was last reset.
 *
 * time: Timer object representing the start time.
 */
long long int ProgressTracker::PassedTime(Timer *tim) {
    chrono::time_point<chrono::high_resolution_clock> clk = chrono::system_clock::now();
    chrono::milliseconds ms = chrono::duration_cast<chrono::milliseconds>(clk - tim->GetClock());

    return ((long long int)ms.count());
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
 *
 * msec:  Number of milliseconds representing the time period.
 * color: Whether or not to print with colors and fancy formatting.
 * lunit: Largest time unit to output.
 */
string ProgressTracker::TimeString(const long long int msec, const bool color, const TimeUnit tu) {
    long long int days, hours, mins, secs, msecs;

    const long long int
        SECOND = 1000,
        MINUTE = 60*SECOND,
        HOUR   = 60*MINUTE,
        DAY    = 24*HOUR;

    if (tu <= DAYS) days = msec/DAY;
    else days = 0;

    if (tu <= HOURS) hours = (msec - days*DAY)/HOUR;
    else hours = 0;

    if (tu <= MINUTES) mins = (msec - days*DAY - hours*HOUR)/MINUTE;
    else mins = 0;

    if (tu <= SECONDS) secs = (msec - days*DAY - hours*HOUR - mins*MINUTE)/SECOND;
    else secs = 0;

    if (tu < MILLISECONDS) msecs = (msec - days*DAY - hours*HOUR - mins*MINUTE - secs*SECOND);
    else msecs = msec;

    return FormattedTimeString(
        days, hours, mins, secs, msecs, color
    );
}

string ProgressTracker::FormattedTimeString(
    const long long int days, const long long int hours, const long long int minutes,
    const long long int seconds, const long long int milliseconds, const bool color
) {
    ostringstream os;
    bool force = false;
    const char *formatstring = "\x1B[1;33m";

    if (days > 0) {
        os << days << "d";
        force = true;
    }

    if (force || hours > 0) {
        if (color) os << formatstring;

        if (hours < 10 && force)
            os << "0";

        if (color)
            os << hours << "\x1B[0mh";
        else
            os << hours << "h";

        force = true;
    }
    
    // MINUTES
    if (force || minutes > 0) {
        if (color) os << formatstring;

        if (minutes < 10 && force)
            os << "0";

        if (color)
            os << minutes << "\x1B[0mm";
        else
            os << minutes << "m";

        force = true;
    }

    // SECONDS
    if (force || seconds > 0) {
        if (color) os << formatstring;

        if (seconds < 10 && force)
            os << "0";
        
        if (color)
            os << seconds << "\x1B[0ms";
        else
            os << seconds << "s";
    }

    // MILLISECONDS
    if (!force && seconds < 10) {
        if (color) os << formatstring;

        if (seconds > 0) {
            if (milliseconds < 10)
                os << "00";
            else if (milliseconds < 100)
                os << "0";
        }

        if (color)
            os << milliseconds << "\x1B[0mms";
        else
            os << milliseconds << "ms";
    }

    return os.str();
}

