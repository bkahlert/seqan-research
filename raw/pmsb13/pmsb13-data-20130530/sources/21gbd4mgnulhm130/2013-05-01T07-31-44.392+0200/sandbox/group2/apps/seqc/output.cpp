#include <iostream>
#include <fstream>

#include "output.h"

Outputter::writeReadLength(const & ReadStats stats)
{    
    // write out stats
    std::ofstream rlf;
    rlf.open("readlength.tsv");
    rlf << "length" << "\t" << "no_of_reads" << std::endl;
    std::map<int,int>::iterator iter; 
   
    for(iter=stats.readLenCount.begin();iter != stats.readLenCount.end(); iter++){
        rlf << iter->first << "\t" << iter->second << std::endl;
    }
    rlf.close();

};

