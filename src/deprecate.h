// Deprecation warning. See: http://stackoverflow.com/questions/295120/c-mark-as-deprecated

#ifndef DEPRECATE_H
#define DEPRECATE_H

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
#pragma message("WARNING: Please implement a DEPRECATED macro for your compiler in MOSGWA file deprecate.h !")
#define DEPRECATED(func) func
#endif

#endif	/* DEPRECATE_H */
