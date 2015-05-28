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
    void writeNucContent(ReadStats & stats, std::stringstream & gc_content,
            std::stringstream & n_content, std::stringstream & nuc_dist);
    
    void writePerReadGcContent(ReadStats & stats, std::stringstream & ss);
    void writePerReadQualContent(ReadStats & stats, std::stringstream & ss);
     
    void writeJobSummary(ReadStats & stats, std::stringstream & ss);
    void guessQualScoreOffset(ReadStats & stats);
};

#endif
