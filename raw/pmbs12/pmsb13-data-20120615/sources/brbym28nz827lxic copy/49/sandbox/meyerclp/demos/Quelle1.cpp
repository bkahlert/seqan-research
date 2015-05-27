#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/find_motif.h>
#include <seqan/file.h>
#include <iostream>

using namespace seqan;



template<typename TAlphabet>
void countAllLetters(TAlphabet const&, String<TAlphabet> str){
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	
	Iterator<String<TAlphabet> >::Type StringIterator = begin(str);
	Iterator<String<TAlphabet> >::Type it2 = end(str);
	while(StringIterator != it2){
		std::cout<<ordValue(*StringIterator)<< " ";
		++StringIterator;
	}


}



int main(){

String<AminoAcid()> str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
//Iterator<String<AminoAcid> >::Type StringIterator = begin(str);
//Iterator<String<AminoAcid> >::Type it2 = end(str);
//while(StringIterator != it2){
//	if(*StringIterator == 'R')
//		*StringIterator = 'A';
//	++StringIterator;
//}
//std::cout<<str;


countAllLetters(AminoAcid(),str);


/*
FrequencyDistribution<AminoAcid> F;
Iterator<String<char> >::Type it1 = begin(str);

absFreqOfLettersInSeq(F,it1,it2);

Iterator<FrequencyDistribution<AminoAcid> >::Type DistributionIterator = begin(F);
Iterator<FrequencyDistribution<AminoAcid> >::Type it4 = end(F);
while(DistributionIterator != it4){
	std::cout<< *DistributionIterator;
	++DistributionIterator;

}
*/


return 0;
}