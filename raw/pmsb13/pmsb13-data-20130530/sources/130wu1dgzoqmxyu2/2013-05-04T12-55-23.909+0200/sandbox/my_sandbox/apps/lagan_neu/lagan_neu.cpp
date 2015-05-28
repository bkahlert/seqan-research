#include <iostream>
#include <fstream>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/align.h>
#include "lagan_functions.h"

using namespace seqan;



int main(int argc, char const ** argv){

	seqan::CharString idOne;
	seqan::DnaString seqOne;
	seqan::CharString idTwo;
	seqan::DnaString seqTwo;

	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq1;
	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq2;

	if (!read_FASTA(argv[1], idOne, seqOne)){
		std::cerr << "error: unable to read first sequence";
		return 1;
	}
	if (!read_FASTA(argv[2], idTwo, seqTwo)){
		std::cerr << "error: unable to read second sequence";
		return 1;
	}
	if (!write_seed_positions(seedsPosSeq1, seedsPosSeq2, argv[3])){
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
	/*typedef Iterator<SSeedSet >::Type SetIterator;

	for (SetIterator it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it){
		std::cout << *it;
		std::cout << std::endl;
	}
	*/


	//typedef Iterator<String<SSeed> > StringIterator;
	seqan::String<SSeed> seedChain;





	unsigned q = 12;
	chainSeedsGlobally(seedChain, seedSet, seqan::SparseChaining());

	for(unsigned i = 0; i < length(seedChain); ++i){
		std::cout << "Initial: seedChain[" << i << "] = " << seedChain[i] << std::endl;
	}

	std::vector<std::pair<unsigned, unsigned> > startPos;
	std::vector<std::pair<unsigned, unsigned> > endPos;

	if (beginPositionH(seedChain[0]) > q && beginPositionV(seedChain[0]) > q){
		startPos.push_back(std::make_pair(0, 0));
		endPos.push_back(std::make_pair(beginPositionH(seedChain[0]), beginPositionV(seedChain[0])));
	}
	for (unsigned i = 0; i < length(seedChain) - 1; ++i){
		unsigned startSeq1 = endPositionH(seedChain[i]);
		unsigned startSeq2 = endPositionV(seedChain[i]);
		unsigned endSeq1 = beginPositionH(seedChain[i + 1]);
		unsigned endSeq2 = beginPositionV(seedChain[i + 1]);

		startPos.push_back(std::make_pair(startSeq1, startSeq2));
		endPos.push_back(std::make_pair(endSeq1, endSeq2));
	}
	if (endPositionH(back(seedChain)) < length(seedChain) - q && endPositionV(back(seedChain)) < length(seedChain) - q){
		startPos.push_back(std::make_pair(endPositionH(back(seedChain)), endPositionV(back(seedChain))));
		endPos.push_back(std::make_pair(length(seqOne), length(seqTwo)));
	}
	unsigned seedChainPosFinder = 0;

	for(unsigned i = 0; i < startPos.size(); ++i){
		if ((endPos[i].first - startPos[i].first) > q && (endPos[i].second - startPos[i].second) > q ){

			seqan::String<SSeed> tempChain;
			if (get_global_seed_chain(tempChain, seqOne, seqTwo, startPos[i], endPos[i], q)){
				insert(seedChain, seedChainPosFinder, tempChain);
				for (unsigned j = 0; j < length(tempChain); ++j){
					std::cout << "i = " << i << "; tempChain[" << j << "] = " << tempChain[j] << std::endl;
				}
			}
			else{
				std::cout << "Keine Treffer in Position " << i << std::endl;
			}
				seedChainPosFinder += length(tempChain) + 1;
			}
		else{
			++seedChainPosFinder;
		}
	}
	//erase(seedChain, 8);
	//erase(seedChain, 16);
	for(unsigned i = 0; i < length(seedChain); ++i){
		std::cout << "seedChain[" << i << "] = " << seedChain[i] << std::endl;
	}

	Align<DnaString, ArrayGaps> alignment;
	resize(seqan::rows(alignment), 2);
	Score<int, Simple> scoringFunction(2, -1, -2);
	seqan::assignSource(seqan::row(alignment, 0), seqOne);
	seqan::assignSource(seqan::row(alignment, 1), seqTwo);

	int result = bandedChainAlignment(alignment, seedChain, scoringFunction, 2);


	//std::cout << "Score: " << result << std::endl;
	//std::cout << alignment << std::endl;

	std::ofstream alignmentOutput;
	alignmentOutput.open ("lagan.txt", std::ios::out | std::ios::trunc | std::ios::binary);
	if (alignmentOutput.is_open()){
		alignmentOutput << "LAGAN: global alignment" << std::endl << std::endl;
		alignmentOutput << "sequence 1: " << idOne << std::endl;
		alignmentOutput << "sequence 2: " << idTwo << std::endl << std::endl;
		alignmentOutput << "score: " << result << std::endl << std::endl;
		alignmentOutput << alignment;
		alignmentOutput.close();
	}
	else{
		std::cerr << "Unable to open file";
	}
	return 0;
}
