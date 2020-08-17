/**
 * Implementation of a simple timing class.
 */

#include <chrono>
#include <string>
#include <sstream>
#include <softlib/Timer.h>

#ifdef UNIX
#   include <sys/time.h>
#endif

using namespace std;
using namespace std::chrono;

/**
 * Constructor.
 */
Timer::Timer() {
    Reset();
}

/**
 * Returns the total elapsed time in milliseconds
 * since this timer was started.
 */
slibreal_t Timer::GetMilliseconds() const {
    struct time t = GetTimeStruct();
    return GetMilliseconds(t);
}
slibreal_t Timer::GetMilliseconds(struct time &t) {
    slibreal_t milliseconds =
        (slibreal_t)(t.tp.count()/1000);

    return milliseconds;
}

/**
 * Returns the total elapsed time in microseconds
 * since this timer was started.
 */
slibreal_t Timer::GetMicroseconds() const {
    struct time t = GetTimeStruct();
    return GetMicroseconds(t);
}
slibreal_t Timer::GetMicroseconds(struct time &t) {
    slibreal_t microseconds =
        (slibreal_t)(t.tp.count());

    return microseconds;
}

/**
 * Returns a Timer::time struct representing the
 * total elapsed time, grouped into days, hours,
 * minutes, seconds, milliseconds and microseconds.
 */
struct Timer::time Timer::GetTimeStruct() const {
    struct time t;

    time_point<high_resolution_clock> toc = high_resolution_clock::now();

    //long long useconds = duration_cast<microseconds>(toc-tic);
    //unsigned long long seconds  = useconds / 1e6;

    t.tp = duration_cast<microseconds>(toc-tic);

    t.days         = duration_cast<hours>(toc-tic) / 24;
    t.hours        = duration_cast<hours>((toc-t.days)-tic);
    t.minutes      = duration_cast<minutes>((toc-t.days-t.hours)-tic);
    t.seconds      = duration_cast<seconds>((toc-t.days-t.hours-t.minutes)-tic);
    t.milliseconds = duration_cast<milliseconds>((toc-t.days-t.hours-t.minutes-t.seconds)-tic);
    t.microseconds = duration_cast<microseconds>((toc-t.days-t.hours-t.minutes-t.seconds-t.milliseconds)-tic);

    return t;
}

/**
 * Reset the timer.
 */
void Timer::Reset() {
    tic = high_resolution_clock::now();
}

/**
 * Returns the currently elapsed time as a string.
 */
string Timer::ToString() const {
    struct time t = GetTimeStruct();
    return ToString(t);
}
string Timer::ToString(struct time &t) {
    ostringstream oss;

    if (t.days.count() > 0)
        oss << t.days.count() << " days, ";
    if (t.hours.count() > 0)
        oss << t.hours.count() << " hours, ";
    if (t.minutes.count() > 0)
        oss << t.minutes.count() << " minutes, ";
    if (t.seconds.count() > 0)
        oss << t.seconds.count() << " seconds, ";
    if (t.milliseconds.count() > 0)
        oss << t.milliseconds.count() << " milliseconds, ";

    oss << t.microseconds.count() << " microseconds";

    return oss.str();
}
