#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>

#include <iostream>
#include <fstream>
#include <map>
#include "output.h"
#include "read_stats.h"


int TSVWriter::writeAll(ReadStats & stats){
    if (mkdir("testout",0744) != 0) {
         std::cerr <<  "Error creating output directory: " << strerror( errno ) << std::endl;
        
        exit(-1);
    } 

    ofStream.open("testout/readlength.tsv");
    
    writeReadLength(stats);    
    ofStream.close();
    
}

int TSVWriter::writeReadLength(ReadStats & stats)
{    
    // write out stats
    ofStream << ">>Sequence Length Distribution\tpass";
    ofStream << "#Length\tcount" << std::endl;
    std::map<int,int>::iterator iter; 
  
    std::cout <<  stats.readLenCount.size();
    for(iter=stats.readLenCount.begin();iter != stats.readLenCount.end(); iter++){
        ofStream << iter->first << "\t" << iter->second << std::endl;
    }
    
    

};


