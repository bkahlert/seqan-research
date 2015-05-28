#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;


int main()
{
    // Build strings
    DnaString str0 = "TATA";
    DnaString str1 = "CGCG"; 
    DnaString str2 = "TTAAGGCC"; 
    DnaString str3 = "ATGC"; 
    DnaString str4 = "AGTGTCA"; 

    // Your code
    StringSet<DnaString> stringSet;
 
    appendValue(stringSet, str0);
    appendValue(stringSet, str1);
    appendValue(stringSet, str2);
    appendValue(stringSet, str3);
     
    
    for(unsigned i = 0; i < length(stringSet); ++i){
    	for(unsigned j = 0; j < length(stringSet[i]); ++j){ std::cout << stringSet[i][j]; }
    	std::cout << std::endl;
    }

    return 0;
}