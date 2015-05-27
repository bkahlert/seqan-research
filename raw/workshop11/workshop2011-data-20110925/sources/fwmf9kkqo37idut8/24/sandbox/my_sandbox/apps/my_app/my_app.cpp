#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout
#include <seqan/basic.h>
using namespace seqan;
using namespace std;

void tut01();
void tut02();
void ass01();

int main(int, char **) {
	ass01();

	return 0;
}

template <typename TAlphabet>
void showAllLetterOfMyAlphabet(TAlphabet const &)
{
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	for (TSize i = 0; i < alphSize; ++i)
		std::cout << i << ',' << TAlphabet(i) << "  ";
	std::cout << std::endl;
}

void tut01(){
	showAllLetterOfMyAlphabet(AminoAcid());
	showAllLetterOfMyAlphabet(Dna());
	showAllLetterOfMyAlphabet(Dna5());
}

void tut02(){
	String<char> str = "ACME"; 
	Iterator<String<char> >::Type it1 = begin(str); //a standard iterator 
	Iterator<String<char>, Standard>::Type it2 = begin(str);  //same as above
	Iterator<String<char>, Rooted>::Type it3 = begin(str);  //a rooted iterator 
	Iterator<String<char>, Rooted>::Type it4 = begin(str, Rooted());  //same as above
}

void ass01(){
	String<AminoAcid> pep = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
	Iterator<String<AminoAcid>, Standard>::Type it = begin(pep);
	Iterator<String<AminoAcid>, Standard>::Type itEnd = end(pep);

	cout << pep;
	while (it != itEnd)
		if (*it == 'R') *it = 'A';
	cout << pep;