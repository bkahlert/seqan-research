#ifndef __OUTPUT_H
#define __OUTPUT_H
#include "read_stats.h"

struct TSVWriter{
   
    const char[] SEQC_END_MODULE_LN = ">>END MODULE"; 
    int writeReadLength(ReadStats & stats);
    int writeAll(ReadStats & stats, char * outdir);

    private:
        std::ofstream ofStream;

        
};

#endif
