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
    time_t start;
    start = time (NULL);
    for (unsigned i = 0; i <= count; ++i)
    {
        appendValue(exact,'A',Exact());
    }
    std::cout << "Strategy Exact() took: " << time(NULL) - start << " s\n\n";
    start = time(NULL);
    for (unsigned i = 0; i <= count; ++i)
        {
            appendValue(generous,'A',Generous());
        }
    std::cout << "Strategy Generous() took: " << time(NULL) - start << " s\n\n";
    return 0;
}
