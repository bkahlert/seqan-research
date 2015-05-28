#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<char> text = "This is the first example";
    Index<String<char>, FMIndex<> > index(text);

    return 0;
}
