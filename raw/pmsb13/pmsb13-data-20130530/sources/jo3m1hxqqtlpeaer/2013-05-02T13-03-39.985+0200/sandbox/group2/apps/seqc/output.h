#ifndef __OUTPUT_H
#define __OUTPUT_H
#include <sstream>
#include "read_stats.h"

struct TSVWriter{
  
    int writeAll(ReadStats & stats, std::string const outdir);
   
    std::stringstream & writeReadLength(ReadStats & stats);
    std::stringstream & writePerPosQual(ReadStats & stats);

    private:
        std::ofstream ofStream;

        
};

#endif
