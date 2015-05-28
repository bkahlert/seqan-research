#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <time.h>

using namespace seqan;

int main()
{
    unsigned count = 1000000;
    String<char> exact;
    String<char> generous;
    clock_t begin=clock();
    begin = clock();
    for (unsigned i = 0; i <= count; ++i)
    {
        appendValue(exact,'A',Exact());
    }
    std::cout << "Strategy Exact() took: " << begin << " s\n\n";
    begin = clock();
    for (unsigned i = 0; i <= count; ++i)
        {
            appendValue(generous,'A',Generous());
        }
    std::cout << "Strategy Generous() took: " << begin << " s\n\n";
    return 0;
}
