#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	
    String<Dna5> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
    Index<String<Dna5>, IndexEsa<> > index(genome);
	Finder<Index<String<Dna5>, IndexEsa<> >> esaFinder(index);
	while(find(esaFinder,"TATAA")){
		std::cout << position(esaFinder) << std::endl;}


	
	StringSet<String<Dna> > dna2;
	appendValue(dna2,  "AC");
	appendValue(dna2,  "CT");
	appendValue(dna2,  "GA");

	Index<StringSet<String<Dna> >, IndexEsa<>> fmIndex(dna2); 

	/*for(unsigned i= 0; i<legth(dna2),++i){
		
	}*/

	
	typedef Index<CharString> TIndex;
    TIndex index2("tobeornottobe");
    Iterator< TIndex, TopDown<ParentLinks<Preorder> > >::Type it(index2);
    while(!atEnd(it))
    {
        std::cout << representative(it) << std::endl;
        ++it;
	}
    return 0;
}