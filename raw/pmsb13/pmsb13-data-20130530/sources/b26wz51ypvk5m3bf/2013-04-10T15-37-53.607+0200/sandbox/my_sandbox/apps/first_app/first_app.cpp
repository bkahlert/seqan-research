#include <seqan/seq_io.h>

int main(int argc, char * argv){
int begin,end;
char file;
char chr;
if(argc==5){
file=argv[1];
chr=argv[2];
begin<<std::istringstream(argv[3]);
end<<std::istringstream(argv[4]);
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
