#include <seqan/sequence.h>
#include <seqan/index.h>
#include <iostream>

using namespace std;
using namespace seqan;

int main()
{
    String<char> text = "This is the first example";
    String<char> pattern = "th";
    Index<String<char>, IndexEsa< > > esaIndex(text);
    Finder<Index<String<char>, IndexEsa<> > > esaFinder(esaIndex);

    find(esaFinder, pattern);
    cout << position(esaFinder) << endl;

    return 0;
}
