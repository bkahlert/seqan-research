#include <seqan/sequence.h>
#include <seqan/index.h>
#include <iostream>

using namespace std;
using namespace seqan;

int main()
{
    String<char> text = "This is the first example";
    Index<String<char>, IndexEsa<> > index(text);

    String<Dna> genome= "TTATTAAGCGTATAGCCCTATAAATATAA";
    Index<String<Dna>, IndexEsa<> > esaIndex(genome); 
    Finder<Index<String<Dna>, IndexEsa<> > > esaFinder(esaIndex);
    
   
   while (find(esaFinder, "TATAA")){
     cout << position(esaFinder) << endl;
   }
   
   CharString texttext="tobeornottobe";
   Index<CharString, IndexEsa<> > esaFinder2(texttext);
   
   Iterator<Index<CharString>, TopDown<ParentLinks<> >  >::Type it(esaFinder2);
   

    do {
      do {
	do{
	  cout << representative(it) << endl;
	} while (goDown(it));
      }while (goRight(it));
    }while (goUp(it)&&!goRight(it));
   
    return 0;
}