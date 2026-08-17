#ifndef GBWTGRAPH_ERROR_HANDLING_H
#define GBWTGRAPH_ERROR_HANDLING_H

/*
  error_handling.h: GBWTGRAPH_THROW, the macro GBWTGraph uses to report
  unrecoverable errors.
*/

#include <stdexcept>

#define GBWTGRAPH_THROW(exception) throw (exception)

#endif // GBWTGRAPH_ERROR_HANDLING_H
