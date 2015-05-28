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

bool get_global_seed_chain(seqan::String<Seed<Simple> > &seedChain, seqan::DnaString &seq1, seqan::DnaString &seq2, unsigned q){

	typedef seqan::SeedSet<Simple> SSeedSet;
	typedef seqan::Seed<Simple> SSeed;
	SSeedSet set;

	seqan::Score<int, seqan::Simple> scoringScheme(1, -1, -1);

	typedef Index<DnaString, IndexQGram<SimpleShape > > qGramIndex;
	qGramIndex index(seq1);

	resize(indexShape(index), q);

	Finder<qGramIndex> myFinder(index);
	std::cout << length(seq2) << std::endl;

	for (unsigned i = 0; i < length(seq2) - (q - 1); ++i){
		DnaString qGram = infix(seq2, i, i + q);
		std::cout << "qGram: " << qGram << std::endl;

		while (find(myFinder, qGram)){
			std::cout << position(myFinder) << std::endl;
			SSeed seed(position(myFinder), i, position(myFinder) + q, i + q);
			if (!addSeed(set, seed, 2, 2, scoringScheme, seq1, seq2, seqan::Chaos())){
				addSeed(set, seed, seqan::Single());
			}
		}
		clear(myFinder);
	}

	if (length(set) > 0){
		seqan::chainSeedsGlobally(seedChain, set, seqan::SparseChaining());
		return true;
	}
	else{
		return false;
	}
}

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

	//creation of seeds and adding them to a SeedSet
	for (unsigned i = 0; i < seedsPosSeq1.size(); ++i){
		SSeed seed(seedsPosSeq1[i].first, seedsPosSeq2[i].first, seedsPosSeq1[i].second, seedsPosSeq2[i].second);
		addSeed(seedSet, seed, Single());
	}
	//std::cout << "Trennlinie" << std::endl;
	typedef Iterator<SSeedSet >::Type SetIterator;

	for (SetIterator it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it){
		std::cout << *it;
		std::cout << std::endl;
	}



	typedef Iterator<String<SSeed> > StringIterator;
	seqan::String<SSeed> initSeedChain;
	std::vector<std::pair<unsigned, seqan::String<SSeed> > > recursionContainer;

	chainSeedsGlobally(initSeedChain, seedSet, seqan::SparseChaining());

	const unsigned minQGram = 12;
	const unsigned intervall = 2;
	const unsigned initQGram = 40 + intervall;
	recursionContainer.push_back(std::make_pair(initQGram, initSeedChain));

	seqan::String<SSeed> globalSeedChain;
	seqan::append(globalSeedChain, initSeedChain[0]);

	while (!recursionContainer.empty()){

		unsigned q = recursionContainer.back().first - intervall;
		int startSeq1 = seqan::endPositionH(seqan::back(globalSeedChain));
		int startSeq2 = seqan::endPositionV(seqan::back(globalSeedChain));

		while (q >= minQGram){

			int endSeq1 = seqan::beginPositionH(seqan::front(recursionContainer.back().second));
			int endSeq2 = seqan::beginPositionV(seqan::front(recursionContainer.back().second));

			if ((endSeq1 - startSeq1) > 5 * q || (endSeq2 - startSeq2) > 5 * q){

				seqan::String<SSeed> tempChain;

				seqan::DnaString window1 = seqan::infix(seqOne, startSeq1, endSeq1);
				seqan::DnaString window2 = seqan::infix(seqTwo, startSeq2, endSeq2);

				std::cout << "window1: " << window1 << std::endl;
				std::cout << "window2: " << window2 << std::endl;

				if (get_global_seed_chain(tempChain, window1, window2, q)){
					std::cout << "get_global_seed_chain: " << length(tempChain) <<  " Treffer!";
					recursionContainer.push_back(std::make_pair(q, tempChain));
				}
				else{
					std::cout << "get_global_seed_chain: Keine Treffer!";
				}
				//last element on a q level
				if (q == minQGram && length(recursionContainer.back().second) == 1){
					append(globalSeedChain, recursionContainer.back().second[0]);
					recursionContainer.pop_back();
					break;
				}
				//no hits on smallest q
				if (q == minQGram && q != recursionContainer.back().first){
					append(globalSeedChain, recursionContainer.back().second[0]);
					erase(recursionContainer.back().second, 0);
					break;
				}
				//hits on smallest q
				if (q == minQGram && recursionContainer.back().first == minQGram){
					seqan::append(globalSeedChain, recursionContainer.back().second);
					recursionContainer.pop_back();
					seqan::append(globalSeedChain, recursionContainer.back().second[0]);
					seqan::erase(recursionContainer.back().second, 0);
					break;
				}
				q -= intervall;
			}
			q -= intervall;
		}

	}
		/*
		std::cout << front(seedChain) << std::endl;
		seqan::append(seedChain2, seedChain[0]);
		std::cout << "bla0" << std::endl;
		seqan::append(seedChain2, seedChain[1]);
		std::cout << length(seedChain2) << std::endl;
		seqan::insert(seedChain2, 1, seedChain);
		std::cout << "seedChain[3]: " << seedChain[3] << std::endl;
		std::cout << "seedChain2[1]: " << seedChain2[1] << std::endl;
		std::cout << "length(seedChain2): " << length(seedChain2) << std::endl;
		erase(seedChain2, 0);
		std::cout << length(seedChain2) << std::endl;
		std::cout << "bla2" << std::endl;

		}
	}

	std::cout << "test2" << std::endl;
	Align<Dna5String, ArrayGaps> alignment;
	resize(seqan::rows(alignment), 2);
	Score<int, Simple> scoringFunction(2, -1, -2);
	seqan::assignSource(seqan::row(alignment, 0), seqOne);
	seqan::assignSource(seqan::row(alignment, 1), seqTwo);



	int result = bandedChainAlignment(alignment, seedChain, scoringFunction, 2);


	std::cout << "Score: " << result << std::endl;
	std::cout << alignment << std::endl;
	*/
		return 0;
}
