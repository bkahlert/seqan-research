#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<char> text = "This is the first example";
    String<char> pattern = "th";
    Index<String<char>, IndexEsa< > > esaIndex(text);
    Finder<Index<String<char>, IndexEsa<> > > esaFinder(esaIndex);

    find(esaFinder, pattern);
    position(esaFinder);

    return 0;
}
