#ifndef GBWTGRAPH_ERROR_HANDLING_H
#define GBWTGRAPH_ERROR_HANDLING_H

/*
  error_handling.h: GBWTGRAPH_THROW, the macro GBWTGraph uses to report
  unrecoverable errors, so that the error-reporting mechanism can be selected
  at build time.

  By default, GBWTGRAPH_THROW(exception_type, message) throws
  exception_type(message), matching GBWTGraph's traditional behavior. Some
  environments (notably Bazel builds embedding GBWTGraph in exceptions-free
  code, such as Google's DeepVariant) need to build with -fno-exceptions.
  Defining GBWTGRAPH_NO_EXCEPTIONS at build time switches GBWTGRAPH_THROW to
  a non-throwing fatal-error path instead: it logs the message with Abseil's
  ABSL_LOG(FATAL) if GBWTGRAPH_USE_ABSEIL_LOGGING is also defined, or
  otherwise prints it to stderr and calls std::abort(). Either way,
  exception_type is not evaluated, so it does not need to be a complete type
  when exceptions are disabled.
*/

#if defined(GBWTGRAPH_NO_EXCEPTIONS)

#if defined(GBWTGRAPH_USE_ABSEIL_LOGGING)
#include "absl/log/absl_log.h"
#define GBWTGRAPH_THROW(exception_type, message) ABSL_LOG(FATAL) << (message)
#else
#include <cstdlib>
#include <iostream>
#define GBWTGRAPH_THROW(exception_type, message) \
    do { std::cerr << (message) << std::endl; std::abort(); } while (0)
#endif

#else

#include <stdexcept>
#define GBWTGRAPH_THROW(exception_type, message) throw exception_type(message)

#endif

#endif // GBWTGRAPH_ERROR_HANDLING_H
