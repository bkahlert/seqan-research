#include <iostream>
#include <seqan/align.h>

using namespace seqan;
int main()
{
	Align< String<AminoAcid> > ali;
    appendValue(rows(ali), "PNCFDAKQRTASRPL");
    appendValue(rows(ali), "CFDKQKNNRTATRDTA");

	LocalAlignmentFinder<> finder(ali);
    Score<int> scoring(3, -2, -1, -5);
	unsigned count = 0;
    while (localAlignment(ali, finder, scoring, 0, WatermanEggert()) && count < 3) {
        ::std::cout << "Score = " << getScore(finder) << ::std::endl;
        ::std::cout << ali;
        ::std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali, 0)) << ":" << (clippedEndPosition(row(ali, 0))-1) << "]";
        ::std::cout << " and Seq2[" << clippedBeginPosition(row(ali, 1)) << ":" <<  (clippedEndPosition(row(ali, 1))-1) << "]" << ::std::endl << ::std::endl;
		++count;
	}
	return 0;
}