#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>
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
    writeReadLength(stats);    
    ofStream.close();
  


}

std::stringstream TSVWriter::writeSomething(ReadStats & stats){
    std::stringstream ss("hallo");
    return ss;
};


int TSVWriter::writeReadLength(ReadStats & stats)
{    
    // write out stats
    ofStream << ">>Sequence Length Distribution\tpass" << std::endl;
    ofStream << "#Length\tCount" << std::endl;
    std::map<int,int>::iterator iter; 
  
    std::cout <<  stats.readLenCount.size();
    for(iter=stats.readLenCount.begin();iter != stats.readLenCount.end(); iter++){
        ofStream << iter->first << "\t" << iter->second << std::endl;
    }
    
    ofStream << SEQC_END_MODULE_LN << std::endl; 

};

int TSVWriter::writePerPosQual(ReadStats & stats)
{
    //>>Per base sequence quality	fail
    //#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile

    //"
};
