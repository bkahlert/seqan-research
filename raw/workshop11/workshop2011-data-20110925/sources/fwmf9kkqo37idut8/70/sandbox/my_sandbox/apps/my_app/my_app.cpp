#include <iostream>
#include <vector>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout
#include <seqan/basic.h>
#include <seqan/find_motif.h>
using namespace seqan;
using namespace std;

void tut01();
void tut02();
void ass01();
void task1();
void countOneMers(String<char> str);

int main(int, char **) {
	task1();

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

	int alphSize = ValueSize<AminoAcid>::VALUE;	
	cout << alphSize << " " << AminoAcid('Q') << endl;

	cout << pep << endl;
	while (it != itEnd){
		
		if (*it == 'R'){
			*it = 'A';
		}
		it++;
	}
	cout << pep << endl;
}

void task1(){
	String<char> str1 = "hello world";
	String<char> str2 = "banana";
	String<char> str3 = "mississippi";

	countOneMers(str1);
	countOneMers(str2);
	countOneMers(str3);
}

template <typename T>
void countOneMers(String<T> str){
	int alphSize = ValueSize<T>::VALUE;
	vector<int> dist (alphSize, 0);

	for (int i=0; i<length(str); i++){
		dist[ordValue(str[i])]++;
	}

	cout << str << endl;
	for (int i=0; i<dist.size(); i++){
		if (dist[i] != 0)
			cout << T(i) << ":" << dist[i] << ", ";
	}
	cout << endl;
}