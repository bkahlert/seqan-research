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
	std::fstream fasta(argv[1], std:: ios_base::in | std::ios_base::binary);
	if(!fasta.good())
		return 1;
	//RecordReader<std::fstream, SinglePass<> > reader(fasta);


	//Variable um Sequenz und ID zu speichern
	typedef String<char,MMap<> > TMMapString;
	TMMapString mmapString;
	if(!open(mmapString, argv[1], OPEN_RDONLY))
		return 1;
	RecordReader<TMMapString, DoublePass<Mapped> > reader(mmapString);


	AutoSeqStreamFormat formatTag;
	if(!checkStreamFormat(reader,formatTag)){
		std::cerr<<"Could not determine file format!"<<std::endl;
		return 1;
	}
	std::cout<<"File format is "<<getAutoSeqStreamFormatName(formatTag)<<'\n';


	StringSet<CharString> ids;
	StringSet<String<Dna5Q> > seqs;

	if(read2(ids,seqs,reader,formatTag) !=0){
 				std::cerr<<"ERROR reading FASTA"<<std::endl;
 				return 1;
 	}
			
	typedef Iterator<StringSet<CharString>, Rooted>::Type TIdIter;
	typedef Iterator<StringSet<String<Dna5Q> >, Standard>::Type TSeqIter;
	TIdIter idIt =begin(ids, Rooted());
	TSeqIter seqIt=begin(seqs, Standard());

	for(;!atEnd(idIt);++idIt,++seqIt){

		
		std::cout<<*idIt<<'\t'<<*seqIt<<std::endl;

	}





	Index<StringSet<String<Dna5Q> > >  index(seqs);
	String<Dna5Q> needle = "A";
	Finder<Index<StringSet<String<Dna5Q> > > > finder(seqs);

	String<char> beg;
	String<char> end;


	Pattern<String<Dna5Q> > pattern(needle);
	while(find(finder, pattern)){
		std::cout<<'[' <<beginPosition(finder)<<','<<endPosition(finder)<<")\t"<<infix(finder)<<std::endl;
		append(beg,beginPosition(finder));
		append(end,endPosition(finder));
	}

	std::cout<<beg<<std::endl;
	std::cout<<end<<std::endl;

	return 0;

}