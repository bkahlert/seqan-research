#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include "output.h"
#include "read_stats.h"

static const char * SEQC_END_MODULE_LN = ">>END_MODULE"; 
static const char * SEQC_READLENGTH_FILE = "readlength.tsv";

 


int TSVWriter::writeAll(ReadStats & stats, std::string const outdir){
    if (mkdir("testout",0744) != 0) {
         std::cerr <<  "Error creating output directory: " << strerror( errno ) << std::endl;
        
        exit(-1);
    } 

    ofStream.open((outdir + "/" + SEQC_READLENGTH_FILE).c_str());
    std::stringstream ss;
    writeReadLength(stats, ss);    
    ofStream << ss;
    ofStream.close();
  


}


void TSVWriter::writeReadLength(ReadStats & stats, std::stringstream & ss)
{    
    // write out stats
    ss << ">>Sequence Length Distribution\tpass" << std::endl;
    ss << "#Length\tCount" << std::endl;
    std::map<int,int>::iterator iter; 
  
    for(iter=stats.readLenCount.begin();iter != stats.readLenCount.end(); iter++){
        ss << iter->first << "\t" << iter->second << std::endl;
    }
    
    ss << SEQC_END_MODULE_LN << std::endl; 

};

void TSVWriter::writePerPosQual(ReadStats & stats, std::stringstream & ss)
{
    //>>Per base sequence quality	fail
    //#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile

    //"
};
