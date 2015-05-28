#include <seqan/basic.h>
#include <seqan/file.h>
#include <impl.h>
#include <seqan/sequence.h>
using namespace seqan;
double quadrat(double d)
{
   return d*d;
}

void iota (String<int> & string, int begin, int end){
    clear(string);
    for(int i=begin;i<end;i++){
    appendValue(string,i);
    }
}

