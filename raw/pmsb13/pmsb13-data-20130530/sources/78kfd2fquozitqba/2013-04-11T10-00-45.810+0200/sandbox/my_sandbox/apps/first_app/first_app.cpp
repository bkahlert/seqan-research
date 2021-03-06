#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/align.h>
#include <time.h>

using namespace seqan;

double diffclock(clock_t clock1, clock_t clock2)
{
    double diffticks = clock1 - clock2;
    double diffms = (diffticks * 1000) / CLOCKS_PER_SEC;
    return diffms;
}
int main()
{
    clock_t begin;
    clock_t end;
    begin = clock();
    String<Dna> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
    String<Dna> pattern1 = "TATAA";
    String<Dna> pattern2 = "TATAA";
    String<Dna> pattern3 = "TATAA";
    StringSet<String<Dna> > pattern;
    appendValue(pattern,pattern1);
    appendValue(pattern,pattern2);
    appendValue(pattern,pattern3);
    Index<StringSet<Dna>, IndexEsa<> > index(text);
    Finder<Index<StringSet<Dna>, IndexEsa<> > > finder(text);
    typedef Align<String<Dna>,ArrayGaps> TAlign;
    TAlign align;
    resize(rows(align),2);
    while (find(finder, pattern))
    {
        std::cout << position(finder) << " " << std::flush;
        Infix<String<Dna> >::Type inf = infix(
                text, position(finder), position(finder) + length(pattern));
//        assignSource(row(align,0),inf);
//        assignSource(row(align,1),pattern);
//        std::cout << align << std::endl;
    }
    end = clock();
    std::cout << "Strategy Generous() took: " << double(diffclock(end, begin))
              << " ms\n\n";
    return 0;
}
