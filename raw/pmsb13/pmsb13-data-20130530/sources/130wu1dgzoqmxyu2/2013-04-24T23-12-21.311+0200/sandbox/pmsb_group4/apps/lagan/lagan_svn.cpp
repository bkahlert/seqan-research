#include <iostream>
#include <fstream>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/seeds.h>

using namespace std;
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

string get_file_contents(const char* filepath){

	ifstream in(filepath, ios::in | ios::binary);
	if (in){
		string contents;

		in.seekg(0, ios::end);
		int fileLength = in.tellg();
		contents.resize(fileLength);
		in.seekg(0, ios::beg);

		in.read(&contents[0], contents.size());
		in.close();
		return(contents);
	}
	throw(errno);
}

bool writeSeedPositions(vector< pair<unsigned, unsigned> > &array1, vector< pair<unsigned, unsigned> > &array2, const char* filepath){

	std::ifstream FileTest(filepath);
	if(!FileTest){
		return false;
	}
	FileTest.close();

	string s;
	s = get_file_contents(filepath);
	string pattern1 = ("Database positions: ");
	string pattern2 = ("Query positions: ");
	pair<unsigned, unsigned> startEnd;

	vector<string> patternContainer;
	patternContainer.push_back(pattern1);
	patternContainer.push_back(pattern2);

	for (unsigned j = 0; j < patternContainer.size(); ++j){

		unsigned found = s.find(patternContainer[j]);

		while (found!=string::npos){
			string temp1 = "";
			string temp2 = "";

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

			stringstream as(temp1);
			stringstream bs(temp2);
			int a;
			int b;
			as >> a;
			bs >> b;
			startEnd = make_pair(a, b);

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

	vector< pair<unsigned, unsigned> > seedsPosSeq1;
	vector< pair<unsigned, unsigned> > seedsPosSeq2;

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
	SSeed seed;

	//setBeginPositionH(seed, seedsPosSeq[0].first);


	std::cout << idOne << std::endl << seqOne << std::endl;
	std::cout << idTwo << std::endl << seqTwo << std::endl;
	cout << "Seq1 " << seedsPosSeq1[0].first << '\t' << seedsPosSeq1[0].second << endl;
	cout << "Seq1 " << seedsPosSeq1[1].first << '\t' << seedsPosSeq1[1].second << endl;
	cout << "Seq2 " << seedsPosSeq2[0].first << '\t' << seedsPosSeq2[0].second << endl;
	cout << "Seq2 " << seedsPosSeq2[1].first << '\t' << seedsPosSeq2[1].second;
	return 0;
}
