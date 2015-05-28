#include <seqan/seq_io.h>

int main(int argc, char** argv){
int begin,end;
char* file;
char* chr;
if(argc==5){
file=argv[1];
chr=argv[2];
begin=seqan::lexicalCast2(begin,argv[4]);
end=seqan::lexicalCast2(end,argv[5]);
}

seqan::FaiIndex faiIdx;
int res = build(faiIdx, file);
if (res != 0)
    std::cerr << "ERROR: Could not build the index!\n";
seqan::CharString str;
unsigned idx = 0;
if (!seqan::getIdByName(faiIdx, chr, idx))
    std::cerr << "ERROR: FAI index has no entry for chr1.\n";
if (readRegion(str, faiIdx, idx, begin, end) != 0)
    std::cerr << "ERROR: Could not load chr1.\n";
std::cout<<str<<std::endl;
return 0;
}
