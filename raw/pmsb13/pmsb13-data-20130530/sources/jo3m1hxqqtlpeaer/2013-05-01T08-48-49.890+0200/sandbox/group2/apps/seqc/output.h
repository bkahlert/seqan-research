#ifndef __OUTPUT_H
#define __OUTPUT_H
#include "read_stats.h"

struct TSVWriter{
   
    int writeReadLength(ReadStats & stats);
    void writeSomething();
    int writeAll(ReadStats & stats);

};

#endif
