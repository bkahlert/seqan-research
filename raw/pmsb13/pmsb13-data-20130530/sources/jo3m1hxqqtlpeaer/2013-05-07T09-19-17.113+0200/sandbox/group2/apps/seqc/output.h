#ifndef __OUTPUT_H
#define __OUTPUT_H
#include <sstream>
#include "read_stats.h"

struct TSVWriter{
  
    int writeAll(ReadStats & stats, std::string const outdir);
   
    void writeReadLength(ReadStats & stats, std::stringstream & ss);
    void writePerPosQual(ReadStats & stats, std::stringstream & ss);

    /**
     * writes out GC content from given ReadStats info into the given stream
     */
    void writeGcContent(ReadStats & stats, std::stringstream & gc_content,
            std::stringstream n_content);
        
};

#endif
