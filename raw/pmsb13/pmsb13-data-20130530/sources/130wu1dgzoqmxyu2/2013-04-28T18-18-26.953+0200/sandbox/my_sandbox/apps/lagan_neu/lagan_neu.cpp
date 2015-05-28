#include <iostream>
#include <fstream>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/align.h>

using namespace seqan;

bool readFASTA(char const * path, CharString  &id, Dna5String &seq){

	std::fstream in(path, std::ios::binary | std::ios::in);
	RecordReader<std::fstream, SinglePass<> > reader(in);

	if (readRecord(id, seq, reader, Fasta()) == 0){
		return true;
	}
	else{
		return false;
	}
}

std::string get_file_contents(const char* filepath){

	std::ifstream in(filepath, std::ios::in | std::ios::binary);
	if (in){
		std::string contents;

		in.seekg(0, std::ios::end);
		int fileLength = in.tellg();
		contents.resize(fileLength);
		in.seekg(0, std::ios::beg);

		in.read(&contents[0], contents.size());
		in.close();
		return(contents);
	}
	throw(errno);
}

bool writeSeedPositions(std::vector< std::pair<unsigned, unsigned> > &array1, std::vector< std::pair<unsigned, unsigned> > &array2, const char* filepath){

	std::ifstream FileTest(filepath);
	if(!FileTest){
		return false;
	}
	FileTest.close();

	std::string s;
	s = get_file_contents(filepath);
	std::string pattern1 = ("Database positions: ");
	std::string pattern2 = ("Query positions: ");
	std::pair<unsigned, unsigned> startEnd;

	std::vector<std::string> patternContainer;
	patternContainer.push_back(pattern1);
	patternContainer.push_back(pattern2);

	for (unsigned j = 0; j < patternContainer.size(); ++j){

		unsigned found = s.find(patternContainer[j]);

		while (found!=std::string::npos){
			std::string temp1 = "";
			std::string temp2 = "";

			int i = 0;
			while (s[found + patternContainer[j].length() + i] != '.'){
				temp1 += s[found + patternContainer[j].length() + i];
				++i;
			}
			i += 2;
			while(s[found + patternContainer[j].length() + i] != '\n'){
						temp2 += s[found + patternContainer[j].length() + i];
						++i;
			}

			std::stringstream as(temp1);
			std::stringstream bs(temp2);
			int a;
			int b;
			as >> a;
			bs >> b;
			startEnd = std::make_pair(a, b);

			if (j == 0){
				array1.push_back(startEnd);
			}
			else{
				array2.push_back(startEnd);
			}
			++found;
			found = s.find(patternContainer[j], found);
		}
	}
	return true;

}


int main(int argc, char const ** argv){


	CharString idOne;
	Dna5String seqOne;
	CharString idTwo;
	Dna5String seqTwo;

	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq1;
	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq2;

	if (!readFASTA(argv[1], idOne, seqOne)){
		std::cerr << "error: unable to read first sequence";
		return 1;
	}
	if (!readFASTA(argv[2], idTwo, seqTwo)){
		std::cerr << "error: unable to read second sequence";
		return 1;
	}
	if (!writeSeedPositions(seedsPosSeq1, seedsPosSeq2, argv[3])){
		std::cerr << "error: STELLAR output file not found";
		return 1;
	}

	typedef Seed<Simple> SSeed;
	typedef SeedSet<Simple> SSeedSet;

	SSeedSet seedSet;

	for (unsigned i = 0; i < seedsPosSeq1.size(); ++i){
		std::cout << "Seq1 " << seedsPosSeq1[i].first << '\t' << seedsPosSeq1[i].second << std::endl;
		std::cout << "Seq2 " << seedsPosSeq2[i].first << '\t' << seedsPosSeq2[i].second << std::endl;
	}

	for (unsigned i = 0; i < seedsPosSeq1.size(); ++i){
		SSeed seed(seedsPosSeq1[i].first, seedsPosSeq2[i].first, seedsPosSeq1[i].second, seedsPosSeq2[i].second);
		std::cout << seed;
		addSeed(seedSet, seed, Single());
		std::cout << "Front: " << beginPositionH(front(seedSet));
		std::cout << std::endl << "back: " << beginPositionH(back(seedSet));
		std::cout << std::endl;
	}

	/*typedef Iterator<SSeedSet >::Type SetIterator;

	for (SetIterator it = begin(seedSet); it != end(seedSet); ++it){
		std::cout << *it;
		std::cout << std::endl;
	}

	*/

	seqan::String<SSeed> seedChain;

	chainSeedsGlobally(seedChain, seedSet, seqan::SparseChaining());

	DnaString window1 = infix(seqOne, endPositionH(seedChain[0]), beginPositionH(seedChain[1]));
	DnaString window2 = infix(seqTwo, endPositionV(seedChain[0]), beginPositionV(seedChain[1]));

	std::cout << "window1: " << window1 << std::endl;
	std::cout << "window2: " << window2 << std::endl;

	//typedef Index<DnaString, IndexQGram<GappedShape<HardwiredShape<1, 2> > > > index("");



	std::cout << "length(seedChain): " << length(seedChain) << std::endl;
	std::cout << "seedChain[0]: " << seqan::beginPositionH(seedChain[0]) << std::endl;
	std::cout << "seedChain[1]: " << seedChain[1] << std::endl;
	std::cout << "seedChain[2]: " << seedChain[2] << std::endl;

	/*std::cout << "test2" << std::endl;
	Align<Dna5String, ArrayGaps> alignment;
	resize(seqan::rows(alignment), 2);
	Score<int, Simple> scoringFunction(2, -1, -2);
	seqan::assignSource(seqan::row(alignment, 0), seqOne);
	seqan::assignSource(seqan::row(alignment, 1), seqTwo);



	int result = bandedChainAlignment(alignment, seedChain, scoringFunction, 2);


	std::cout << "Score: " << result << std::endl;
	std::cout << alignment << std::endl;
	*/


	/*std::cout << idOne << std::endl << seqOne << std::endl;
	std::cout << idTwo << std::endl << seqTwo << std::endl;
	for (unsigned i = 0; i < seedsPosSeq1.size(); ++i){
		cout << "Seq1 " << seedsPosSeq1[i].first << '\t' << seedsPosSeq1[i].second << endl;
		cout << "Seq2 " << seedsPosSeq2[i].first << '\t' << seedsPosSeq2[i].second << endl;
	}
	*/
		return 0;
}
