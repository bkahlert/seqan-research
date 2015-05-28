#include <seqan/seq_io.h>

int main(int argc, char** argv){
int begin,end;
char* file,chr;
if(argc==5){
file=argv[1];
chr=argv[2];
begin=lexicalCast2(begin,argv[4]);
end=lexicalCast2(end,argv[5]);
}

seqan::FaiIndex faiIndex;
int res = build(faiIndex, file);
if (res != 0)
    std::cerr << "ERROR: Could not build the index!\n";
seqan::CharString str;
if (readRegion(str, faiIdx, idx, begin, end) != 0)
    std::cerr << "ERROR: Could not load chr1.\n";
std::cout<<str<<std::endl;
return 0;
}
