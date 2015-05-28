#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<char> text = "This is the first example";
    String<char> pattern = "th";
    Index<String<char>, FMIndex<> > fmIndex(text);
    Finder<Index<String<char>, FMIndex<> > > fmFinder(fmIndex);

    find(fmFinder, pattern);
    position(fmFinder);

    return 0;
}
