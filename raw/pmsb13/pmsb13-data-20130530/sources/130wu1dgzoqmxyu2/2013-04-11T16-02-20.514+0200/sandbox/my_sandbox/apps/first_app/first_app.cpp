#include <iostream>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <vector>

using namespace seqan;
using namespace std;

void assignmentIndex2(){
	String<char> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
	StringSet<String<char> > set;

	String<char> pattern1 = "TA";
	String<char> pattern2 = "GCG";
	String<char> pattern3 = "GCC";

	appendValue(set, pattern1);
	appendValue(set, pattern2);
	appendValue(set, pattern3);

	Index<String<char>, IndexEsa<> > esaIndex(text);
	Finder<Index<String<char>, IndexEsa<> > > esaFinder(esaIndex);

	typedef Iterator<StringSet<String<char> > >::Type TStringSetIterator;
	for (TStringSetIterator it = begin(set); it != end(set); ++it){
		clear(esaFinder);
		std::cout << value(it) << ": ";
		while(find(esaFinder, getValue(it))){
			std::cout<<position(esaFinder) << " ";
		}
		std::cout << "\n";
	}
}

void assignmentIndexIterators2(){
	String<char> text = "tobeornottobe";
	typedef Index<String<char> > TIndex;
	TIndex index(text);

	Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);
	bool stop = false;
	while(!stop){
		cout << representative(it) << endl;
		if (!goDown(it))
			while(!goRight(it)){
				goUp(it);
				if (isRoot(it)){
					stop = true;
					break;
				}
			}
	}
}

void assignment1qGram(){
	typedef Index<DnaString, IndexQGram<OneGappedShape> > TIndex;
	TIndex index("CATGATTACATA");
	stringToShape(indexShape(index), "1101");
	hash(indexShape(index), "AT-A");
	for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i){
		std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

	}
}

void assignment2qGram(){
	vector<vector<int > > matrix;
	matrix.resize(3);
	for (int i = 0; i < 3; ++i){
		matrix[i].resize(3,0);
	}
	StringSet<DnaString> set;
	appendValue(set, "ACTTTATACCTGTATGAACGTAGCTAGTGTTGGATTGATGACACTGATCGATGGTCTTCTACTTCTGTCATGTA");
	appendValue(set, "CCCTGTCGGTGATGTCGTTAGCTATGCGTAGCGTCGTTGCAGTGTCATGTATTCGGCGTAGCGGACTGATCTCGCTAGTCATGACTGTAGCTAGCGCTGTGGCGGGATT");
	appendValue(set, "AAGGGTTCCGATAGCGAG");


	//typedef Index<DnaString, IndexQGram< UngappedShape<5> > > index(set);
	for (int i = 0; i < 3; ++i){
		for (int j = i + 1; j < 3 ; ++j){
			cout << "Moritz ist faul!!!" << endl;
		}
	}
}

int main(){
	assignment2qGram();
	return 0;
}
