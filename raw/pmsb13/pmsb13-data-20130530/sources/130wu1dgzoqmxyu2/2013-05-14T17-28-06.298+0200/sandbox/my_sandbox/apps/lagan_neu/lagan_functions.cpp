#include "lagan_functions.h"

//using namespace seqan;

bool get_global_seed_chain(seqan::String<seqan::Seed<seqan::Simple> > &seedChain,
							seqan::DnaString &seq1, seqan::DnaString &seq2,
							std::pair<unsigned, unsigned> startPos,
							std::pair<unsigned, unsigned> endPos, unsigned q)
{

	seqan::DnaString window1 = seqan::infix(seq1, startPos.first, endPos.first);
	seqan::DnaString window2 = seqan::infix(seq2, startPos.second, endPos.second);

	typedef seqan::SeedSet<seqan::Simple> SSeedSet;
	typedef seqan::Seed<seqan::Simple> SSeed;
	SSeedSet set;

	seqan::Score<int, seqan::Simple> scoringScheme(2, -1, -1);
	typedef seqan::Index<seqan::DnaString, seqan::IndexQGram<seqan::SimpleShape> > qGramIndex;
	qGramIndex index(window1);
	resize(indexShape(index), q);

	seqan::Finder<qGramIndex> myFinder(index);
	//find all common qGrams and add them to a SeedSet (if Chaos() fails, add them Single())
	for (unsigned i = startPos.second; i < seqan::length(window2) + startPos.second - (q - 1); ++i)
	{
		seqan::DnaString qGram = seqan::infix(seq1, i, i + q);
		while (seqan::find(myFinder, qGram))
		{
			/*debugging:
				std::cout << "In Funktion! qGram: " << qGram << std::endl;
				std::cout << "i: " << i << std::endl;
				std::cout << position(myFinder) << std::endl;
			*/
			SSeed seed(position(myFinder) + startPos.first, i, position(myFinder) + q + startPos.first, i + q);
			if (!seqan::addSeed(set, seed, 2, 2, scoringScheme, seq1, seq2, seqan::Chaos()))
			{
				seqan::addSeed(set, seed, seqan::Single());
			}
		}
		seqan::clear(myFinder);
	}
	if (seqan::length(set) > 0)
	{
		seqan::chainSeedsGlobally(seedChain, set, seqan::SparseChaining());
		return true;
	}
	else
	{
		return false;
	}
}

bool insert_seeds_in_seedChain(seqan::String<seqan::Seed<seqan::Simple> > &seedChain,
										std::vector<std::pair<unsigned, unsigned> > &startPos,
										std::vector<std::pair<unsigned, unsigned> > &endPos,
										seqan::DnaString &seqOne, seqan::DnaString &seqTwo,
										unsigned q)
{

	//check if any container is empty
	if(!startPos.empty() && !endPos.empty() &&
		seqan::length(seedChain) > 0 && seqan::length(seqOne) > 0 && length(seqTwo) > 0)
	{

		unsigned seedChainPosFinder = 0; //determines the position of the next, unsearched gap

		//for all gaps between the stellar matches...
		for(unsigned i = 0; i < startPos.size(); ++i)
		{
			std::cout << "i = " << i << " ;seedChainPosFinder = " << seedChainPosFinder << std::endl;
			if ((endPos[i].first - startPos[i].first) > q &&
				(endPos[i].second - startPos[i].second) > q )
			{

				//insert the found seeds between the corresponding stellar matches
				seqan::String<seqan::Seed<seqan::Simple> > tempChain;
				if (get_global_seed_chain(tempChain, seqOne, seqTwo, startPos[i], endPos[i], q))
				{
					seqan::insert(seedChain, seedChainPosFinder, tempChain);
					for (unsigned j = 0; j < seqan::length(tempChain); ++j)
					{
						std::cout << "i = " << i << "; tempChain[" << j << "] = " << tempChain[j] << std::endl;
					}
				}
				else
				{
					std::cout << "Keine Treffer in Gap-Position " << i << std::endl;
				}
				seedChainPosFinder += seqan::length(tempChain) + 1;

			}
			else{
				++seedChainPosFinder;
			}
		}
		return true;
	}
	else{
		return false;
	}
}

bool create_positionsContainers(std::vector<std::pair<unsigned, unsigned> > &startPos,
									std::vector<std::pair<unsigned, unsigned> > &endPos,
									seqan::String<seqan::Seed<seqan::Simple> > &seedChain,
									unsigned q, unsigned lenSeqOne, unsigned lenSeqTwo)
{
	if (seqan::empty(seedChain))
	{
		std::cerr << "Leere seedChain übergeben..." << std::endl;
		return false;
	}
	std::cerr << "bumsen" << std::endl;
	if (seqan::beginPositionH(seedChain[0]) > 2 * q && seqan::beginPositionV(seedChain[0]) > 2 * q)
	{
		startPos.push_back(std::make_pair(0, 0));
		endPos.push_back(std::make_pair(seqan::beginPositionH(seedChain[0]), seqan::beginPositionV(seedChain[0])));
	}
	std::cerr << "ficken" << std::endl;
	for (unsigned i = 0; i < seqan::length(seedChain) - 1; ++i)
	{
		unsigned startSeq1 = seqan::endPositionH(seedChain[i]);
		unsigned startSeq2 = seqan::endPositionV(seedChain[i]);
		unsigned endSeq1 = seqan::beginPositionH(seedChain[i + 1]);
		unsigned endSeq2 = seqan::beginPositionV(seedChain[i + 1]);

		startPos.push_back(std::make_pair(startSeq1, startSeq2));
		endPos.push_back(std::make_pair(endSeq1, endSeq2));
	}
	std::cerr << "blasen" << std::endl;
	if (seqan::endPositionH(back(seedChain)) < lenSeqOne - (2 * q) &&
		seqan::endPositionV(seqan::back(seedChain)) < lenSeqTwo - (2 * q))
	{
		startPos.push_back(std::make_pair(seqan::endPositionH(seqan::back(seedChain)),
											seqan::endPositionV(seqan::back(seedChain))));
		endPos.push_back(std::make_pair(lenSeqOne, lenSeqTwo));
	}
	std::cerr << "alles" << std::endl;
	return true;
}

bool read_FASTA(seqan::CharString path, seqan::CharString  &id, seqan::DnaString &seq)
{

	std::fstream in(seqan::toCString(path), std::ios::binary | std::ios::in);
	seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);

	if (seqan::readRecord(id, seq, reader, Fasta()) == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

std::string get_file_contents(seqan::CharString filepath)
{

	std::ifstream in(seqan::toCString(filepath), std::ios::in | std::ios::binary);
	if (in)
	{
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


bool write_seed_positions(std::vector< std::pair<unsigned, unsigned> > &array1,
							std::vector< std::pair<unsigned, unsigned> > &array2,
							const char* filepath)
{

	std::ifstream FileTest(filepath);
	if(!FileTest)
	{
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

	for (unsigned j = 0; j < patternContainer.size(); ++j)
	{

		unsigned found = s.find(patternContainer[j]);

		while (found!=std::string::npos)
		{
			std::string temp1 = "";
			std::string temp2 = "";

			int i = 0;
			while (s[found + patternContainer[j].length() + i] != '.')
			{
				temp1 += s[found + patternContainer[j].length() + i];
				++i;
			}
			i += 2;
			while(s[found + patternContainer[j].length() + i] != '\n')
			{
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

			if (j == 0)
			{
				array1.push_back(startEnd);
			}
			else
			{
				array2.push_back(startEnd);
			}
			++found;
			found = s.find(patternContainer[j], found);
		}
	}
	return true;

}

bool output_alignment(seqan::CharString path, seqan::CharString &idOne, seqan::CharString &idTwo,
						seqan::Align<seqan::DnaString, seqan::ArrayGaps> &alignment,
						int score)
{
	std::ofstream alignmentOutput;
	alignmentOutput.open (seqan::toCString(path), std::ios::out | std::ios::trunc | std::ios::binary);
	if (alignmentOutput.is_open())
	{
		alignmentOutput << "LAGAN: global alignment" << std::endl << std::endl;
		alignmentOutput << "sequence 1 ID: " << idOne << std::endl;
		alignmentOutput << "sequence 2 ID: " << idTwo << std::endl << std::endl;
		alignmentOutput << "score: " << score << std::endl << std::endl;
		alignmentOutput << alignment;
		alignmentOutput.close();
		return true;
	}
	else
	{
		return false;
	}
}

