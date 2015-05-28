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

	
	//typedef Index<CharString> TIndex;
    Index<CharString> index2("tobeornottobe");
    //Iterator<Index<CharString>, TopDown<ParentLinks<>> >::Type it(index2);
	Iterator< Index<CharString>, TopDown<ParentLinks<Preorder> > >::Type it(index2);
    while(!atEnd(it))
    {
        std::cout << representative(it) << std::endl;
        ++it;
	}
	/*do{
		//std::cout << representative(it) << std::endl;
		while(goDown(it)){
			std::cout << representative(it) << std::endl;
		}
		goRight(it);
	}while(!isRoot(it));*/
    
	return 0;
}