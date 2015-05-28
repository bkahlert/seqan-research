#ifndef __OUTPUT_H
#define __OUTPUT_H
#include "read_stats.h"

struct TSVWriter{
   
    int writeReadLength(ReadStats & stats);
    int writeAll(ReadStats & stats, std::string const outdir);
    int writePerPosQual(ReadStats & stats);
    void writeSomething(ReadStats & stats, std::stringstream & ss);

    private:
        std::ofstream ofStream;

        
};

#endif
