#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
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
    String<Dna> pattern = "TATAA";
    Index<String<Dna>, IndexEsa<> > index(text);
    Finder<Index<String<Dna>, IndexEsa<> >> finder(text);
    while(find(finder,pattern)){
        std::cout << position(finder);
    }
    end = clock();
    std::cout << "Strategy Generous() took: " << double(diffclock(end, begin))
              << " ms\n\n";
    return 0;
}
