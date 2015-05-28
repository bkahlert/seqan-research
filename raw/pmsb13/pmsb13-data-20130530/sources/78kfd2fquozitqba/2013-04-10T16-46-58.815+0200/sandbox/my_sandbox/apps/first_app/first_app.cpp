#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <time.h>

using namespace seqan;

double diffclock(clock_t clock1,clock_t clock2)
{
  double diffticks=clock1-clock2;
  double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
  return diffms;
}
int main()
{
    unsigned count = 1000000;
    String<char> exact;
    String<char> generous;
//    clock_t begin=clock();
//    begin = clock();
//    for (unsigned i = 0; i <= count; ++i)
//    {
//        appendValue(exact,'A',Exact());
//    }
//    std::cout << "Strategy Exact() took: " << begin << " s\n\n";
    clock_t begin=clock();
    for (unsigned i = 0; i <= count; ++i)
        {
            appendValue(generous,'A',Generous());
        }
    clock_t end=clock();
    std::cout << "Strategy Generous() took: " << double(diffclock(end,begin)) << " ms\n\n";
    begin=clock();
        for (unsigned i = 0; i <= count; i++)
            {
                appendValue(generous,'A',Generous());
            }
        end=clock();
    std::cout << "Strategy Generous() took: " << double(diffclock(end,begin)) << " ms\n\n";
    return 0;
}
