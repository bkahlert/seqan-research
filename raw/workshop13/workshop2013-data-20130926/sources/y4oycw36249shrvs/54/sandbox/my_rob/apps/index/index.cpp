#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna> text = "This is the first example";
    Index<String<Dna>, IndexEsa<> > index(text);
    StringSet<String<Dna> > stringSet;
    
    appendValue(stringSet, "This");
    appendValue(stringSet, "That");
	
    Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(stringSet);
    

    return 0;
}