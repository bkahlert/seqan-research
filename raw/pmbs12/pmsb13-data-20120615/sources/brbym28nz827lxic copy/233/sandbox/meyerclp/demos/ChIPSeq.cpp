#include <fstream>
#include <iostream>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/stream.h>
#include <seqan/index.h>
 
using namespace seqan;



int main(int argc, char const ** argv){

	if(argc !=2){
		std::cerr<<"ERROR: Invalid argument count."<<std:: endl
				 <<"Usage:" <<argv[0]<<"File"<<std::endl;
		return 1;
	}
	std:: ifstream fasta(argv[1], std:: ios_base::in | std::ios_base::binary);
	if(!fasta.good())
		return 1;
	RecordReader<std::ifstream, SinglePass<> > reader(fasta);

	AutoSeqStreamFormat formatTag;
	if(!checkStreamFormat(reader,formatTag)){
		std::cerr<<"Could not determine file format!"<<std::endl;
		return 1;
	}
	std::cout<<"File format is "<<getAutoSeqStreamFormatName(formatTag)<<'\n';


	//Variable um Sequenz und ID zu speichern
	CharString id;
	Dna5String seq;

	while(!atEnd(reader)){

		if(readRecord(id,seq,reader,formatTag) !=0){
			std::cerr<<"ERROR reading FASTA"<<std::endl;
			return 1;
		}
		std::cout<<id<<"\t"<<seq<<"\n";
	}
	
	Index<Dna5String> index(seq);
	Dna5String needle = "A";
	Finder<Index<Dna5String> > finder(index);

	Pattern<Dna5String> pattern(needle);
	while(find(finder, pattern)){
		std::cout<<'[' <<beginPosition(finder)<<','<<endPosition(finder)<<")\t"<<infix(finder)<<std::endl;
	}

	return 0;

}