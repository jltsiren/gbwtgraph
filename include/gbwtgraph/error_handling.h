#ifndef GBWTGRAPH_ERROR_HANDLING_H
#define GBWTGRAPH_ERROR_HANDLING_H

/*
  error_handling.h: GBWTGRAPH_THROW, the macro GBWTGraph uses to report
  unrecoverable errors, so that the error-reporting mechanism can be selected
  at build time.

  GBWTGRAPH_THROW(exception) throws the given exception object, unless
  GBWTGRAPH_NO_EXCEPTIONS is defined, in which case it logs exception.what()
  and calls std::abort() instead (via Abseil if GBWTGRAPH_USE_ABSEIL_LOGGING
  is set).
*/

// GBWTGRAPH_THROW always constructs its argument, even when not thrown.
#include <stdexcept>

#if defined(GBWTGRAPH_NO_EXCEPTIONS)

#if defined(GBWTGRAPH_USE_ABSEIL_LOGGING)
#include "absl/log/absl_log.h"
#define GBWTGRAPH_THROW(exception) ABSL_LOG(FATAL) << (exception).what()
#else
#include <cstdlib>
#include <iostream>
#define GBWTGRAPH_THROW(exception) \
    do { std::cerr << (exception).what() << std::endl; std::abort(); } while (0)
#endif

#else

#define GBWTGRAPH_THROW(exception) throw (exception)

#endif

#endif // GBWTGRAPH_ERROR_HANDLING_H
