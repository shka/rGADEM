/* src/config.h.  Generated from config.h.in by configure.  */
#ifndef CONFIG_H
#define CONFIG_H

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "http://www.rglab.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "rGADEM"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "rGADEM 1.1.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "rgadem"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.1.1"

/* Define to 1 if you have the <dispatch/dispatch.h> header file. */
#define HAVE_DISPATCH_DISPATCH_H 1

/* Define to 1 if you have the OpenMP support. */
/* #undef HAVE_OPENMP */

#ifdef HAVE_DISPATCH_DISPATCH_H
  #include <dispatch/dispatch.h>
  #define DO_APPLY(task, n_times, counter_name) \
    R_CheckUserInterrupt(); \
    dispatch_apply(n_times, \
                   dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),	\
                   ^(size_t counter_name) { task; });
  
#else // not HAVE_DISPATCH_DISPATCH_H
  #ifdef HAVE_OPENMP
    #include <omp.h>
    #define DO_APPLY(task, n_times, counter_name) \
      R_CheckUserInterrupt(); \
      _Pragma("omp parallel for") \
      for (int counter_name = 0; counter_name < n_times; ++counter_name) { \
        task; \
      }
  #else // not HAVE_OPENMP
    #define DO_APPLY(task, n_times, counter_name) \
      for (int counter_name = 0; counter_name < n_times; ++counter_name) { \
        R_CheckUserInterrupt(); \
        task; \
      }
  #endif // HAVE_OPENMP
#endif // HAVE_DISPATCH_DISPATCH_H

#endif
