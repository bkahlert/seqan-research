#include <fstream>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/store.h>

using namespace seqan;

int main()
{
    FragmentStore<> store;

    std::ifstream file("assignment_annotations.gtf", std::ios_base::in | std::ios_base::binary);
    read(file, store, Gtf());
    // Create AnnotationTree iterator
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
    // Move iterator one node down 
#if 0
    while (goDown(it)) ;

    int count = 1;
    while (goRight(it))
        ++count;
    std::cout << "n=" << count << std::endl;
#elif 0
    unsigned count = 0;
    while (goDown(it)) ;
    while (!atEnd(it)){
        ++count;
        // Iterate over all siblings and count
        while (goRight(it))
            ++count;
        std::cout << count << std::endl;
        count = 0;
        // Jump to the next mRNA or gene, go down to its first leaf and count it
        if (!atEnd(it)) {
            goNext(it);
            if (!atEnd(it)) while(goDown(it)) ;
        }
    }
#else
    while (goDown(it)) ;
    std::cout << "type: " <<  getType(it) << std::endl;
    std::cout << "begin position: " <<  getAnnotation(it).beginPos << std::endl;
    std::cout << "end position: " <<  getAnnotation(it).endPos << std::endl;
    std::cout << "id: " << value(it) << std::endl;
    goUp(it);
    std::cout << "parent id: " << value(it) << std::endl;
    std::cout << "parent name: " <<  getName(it) << std::endl;

    goRoot(it);
    unsigned no_gene = 0;
    unsigned no_mrna = 0;
    unsigned no_exon = 0;
    unsigned exlen = 0;
    while (!atEnd(it))
    {
	    if (getType(it) == "gene")
		    ++no_gene;
	    if (getType(it) == "mrna")
		    ++no_mrna;
	    if (getType(it) == "exon") {
		    ++no_exon;
		    exlen += getAnnotation(it).endPos - getAnnotation(it).beginPos;
	    }
	    goNext(it);
    }
#endif
    std::cout << "gene " << no_gene << std::endl;
    std::cout << "mrna " << no_mrna << " " << no_mrna / no_gene << std::endl;
    std::cout << "exon " << no_exon << " " << no_exon / no_mrna << std::endl;
    std::cout << "     " << exlen << " " << exlen / no_exon << std::endl;
 
    return 0;
}