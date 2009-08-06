#include <cstdarg>
#include <cstdio>
#include <cstdlib>

namespace lshkit {

    void panic_intern(const char *fmt, ...)
    {
        va_list args;
        char msg[2048];

        va_start(args, fmt);
        vsnprintf(msg, sizeof(msg), fmt, args);
        va_end(args);
        fputs(msg, stderr);
        exit(-1);
    }

}

