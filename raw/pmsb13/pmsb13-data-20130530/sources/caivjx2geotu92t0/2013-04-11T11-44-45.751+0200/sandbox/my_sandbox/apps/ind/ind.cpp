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

	// Q-gram index

	//typedef Index<DnaString, IndexQGram< OneGappedShape<4>> >  TIndex;
	Shape<OneGappedShape> shape;
    Index<DnaString, IndexQGram< Shape<OneGappedShape>> > index3("CATGATTACATA");
    stringToShape(shape, '1101');
	 hash(indexShape(index3), "AT-A");
    for (unsigned i = 0; i < length(getOccurrences(index3, indexShape(index3))); ++i)
        std::cout << getOccurrences(index3, indexShape(index3))[i] << std::endl;
	return 0;
}