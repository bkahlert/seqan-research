#ifndef __OUTPUT_H
#define __OUTPUT_H
#include "read_stats.h"

struct TSVWriter{
    
    int writeReadLength(ReadStats const & stats);
    int writeAll(ReadStats const & stats);

};

#endif
