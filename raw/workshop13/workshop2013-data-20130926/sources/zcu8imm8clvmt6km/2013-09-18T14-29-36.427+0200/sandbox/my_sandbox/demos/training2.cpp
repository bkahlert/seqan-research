#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
	Align< String<char> > ali;
	resize(rows(ali), 2);
	assignSource(row(ali, 0), "aphilologicaltheorem");
	assignSource(row(ali, 1), "bizarreamphibology");
	std::cout << "Score = " << localAlignment(ali, Score<int>(3, -3, -2, -2)) << std::endl;
	std::cout << ali;
	std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali, 0)) << ":" << (clippedEndPosition(row(ali, 0)) - 1) << "]";
	std::cout << " and Seq2[" << clippedBeginPosition(row(ali, 1)) << ":" << (clippedEndPosition(row(ali, 1)) - 1) << "]" << std::endl << std::endl;
	Align< String<Dna> > ali2;
	resize(rows(ali2), 2);
	assignSource(row(ali2, 0), "ataagcgtctcg");
	assignSource(row(ali2, 1), "tcatagagttgc");
	Score<int> scoring(2, -1, -2, 0);
	LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(scoring, 5);
	while (nextLocalAlignment(ali2, enumerator))
	{
		std::cout << "Score = " << getScore(enumerator) << std::endl;
		std::cout << ali2;
		std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali2, 0)) << ":" << (clippedEndPosition(row(ali2, 0)) - 1) << "]";
		std::cout << " and Seq2[" << clippedBeginPosition(row(ali2, 1)) << ":" << (clippedEndPosition(row(ali2, 1)) - 1) << "]" << std::endl << std::endl;
	}
	return 0;
}