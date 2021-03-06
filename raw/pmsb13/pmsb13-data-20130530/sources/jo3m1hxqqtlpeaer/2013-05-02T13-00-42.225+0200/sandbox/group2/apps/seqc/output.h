#ifndef __OUTPUT_H
#define __OUTPUT_H
#include <sstream>
#include "read_stats.h"

struct TSVWriter{
   
    int writeReadLength(ReadStats & stats);
    int writeAll(ReadStats & stats, std::string const outdir);
    int writePerPosQual(ReadStats & stats);
    std::stringstream & writeSomething(ReadStats & stats, std::stringstream & ss);

    private:
        std::ofstream ofStream;

        
};

#endif
