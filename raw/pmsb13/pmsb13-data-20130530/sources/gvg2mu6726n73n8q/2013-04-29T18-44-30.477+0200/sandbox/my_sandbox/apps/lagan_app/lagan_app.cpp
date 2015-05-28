#include <iostream>
#include <seqan/index.h>
#include <vector>
 
using namespace seqan;
using namespace std;
 
 
int main(){
 
    DnaString window1 = "AAAACACGC";
    DnaString window2 = "GGTCGACCGT";
 
    StringSet<DnaString > stringSet;
    appendValue(stringSet, window1);
    appendValue(stringSet, window2);
 
    typedef Index< StringSet<DnaString>, IndexQGram<UngappedShape<2> > > qGramIndex;
    qGramIndex index(stringSet);
    Finder<qGramIndex> myFinder(index);
 
    while (find(myFinder, "CG")){
        std::cout << position(myFinder) << "  ";
    }
 
    std::cout << std::endl << indexCountsDir(index);
 
    return 0;
}