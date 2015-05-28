#ifndef __OUTPUT_H
#define __OUTPUT_H
#include <sstream>
#include "read_stats.h"

struct TSVWriter{
  
    int writeAll(ReadStats & stats, std::string const outdir);
   
    void writeReadLength(ReadStats & stats, std::stringstream & ss);
    void writePerPosQual(ReadStats & stats, std::stringstream & ss);

    private:

        
};

#endif
