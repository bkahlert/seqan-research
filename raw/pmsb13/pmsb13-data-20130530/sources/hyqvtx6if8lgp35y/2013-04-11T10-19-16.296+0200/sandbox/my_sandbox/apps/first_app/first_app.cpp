#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	String<char> text = "TTATTAAGCGTATAGCCCTATAAATATAA";    
	String<char> pattern1 = "TATAA";
	String<char> pattern2 = "GC";
	String<char> pattern3 = "AA";
	StringSet<String<char> > allPattern;
	appendValue(allPattern, pattern1);
	appendValue(allPattern, pattern2);
	appendValue(allPattern, pattern3); 
	
	Index<String<char> > esaIndex(text);
	Finder<Index<String<char> > > esaFinder(esaIndex);
	
	
	for(unsigned i = 0; i < length(allPattern); ++i)
	{	
		while(find(esaFinder, allPattern[i]))
		{		
			std::cout << "Occurence of Pattern " << allPattern[i] << '\t' <<;			
			std::cout << position(esaFinder) <<std::endl;
		}		
		clear(esaFinder);	
	}
    return 0;
}
