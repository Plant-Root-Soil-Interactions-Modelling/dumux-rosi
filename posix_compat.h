#ifndef POSIX_COMPAT_H
#define POSIX_COMPAT_H

#ifdef _WIN32

#include <direct.h>
#include <cstdlib>
#include <cstdio>
#include <io.h>
#include <chrono>
#include <ctime>

// drand48 replacement
inline double drand48()
{
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}
#define srand48(x) srand(x)

// setenv replacement
inline int setenv(const char *name, const char *value, int overwrite)
{
    if (!overwrite) {
        size_t requiredSize = 0;
        errno_t err = getenv_s(&requiredSize, nullptr, 0, name);
        if (err == 0 && requiredSize > 0) return 0; // already set
    }
    return _putenv_s(name, value);
}

#endif // _WIN32
#endif // POSIX_COMPAT_H

