#include <seqan/seq_io.h>

int main(int argc, char const argv){
int begin,end;
char const file,chr;
if(argc==5){
file=argv[1];
chr=argv[2];
begin=argv[3];
end=argv[4];
}

seqan::FaiIndex faiIndex;
int res = build(faiIndex, "path/to/file.fasta");
if (res != 0)
    std::cerr << "ERROR: Could not build the index!\n";
seqan::CharString str;
if (readRegion(str, faiIdx, idx, begin, end) != 0)
    std::cerr << "ERROR: Could not load chr1.\n";
std::cout<<str<<std::endl;
return 0;
}
