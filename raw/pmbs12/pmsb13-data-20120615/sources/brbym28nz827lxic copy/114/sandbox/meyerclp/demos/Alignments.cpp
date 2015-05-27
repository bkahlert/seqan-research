#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/find_motif.h>
#include <seqan/file.h>
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_align.h>

using namespace seqan;

int main (){

	/*typedef String<Dna> TDnaSeq;
	typedef Align<TDnaSeq, ArrayGaps> TAlign;


	TDnaSeq seq1 = "acgtacgtact";
	TDnaSeq seq2 = "actactacgt";
	TAlign align;
	resize(rows(align),2);
	assignSource(row(align,0),seq1);
	assignSource(row(align,1),seq2);

	::std::cout<<align;

	int score = globalAlignment(align,Score<int>(1,-1,-1,-1), Hirschberg());

	::std::cout<<align;

	::std::cout<<row(align,1);

	typedef Row<TAlign>::Type TRow;
	TRow &row1 = row(align,0);
	TRow &row2 = row(align,1);

	Iterator<TRow>::Type AlignIter = end(row1);
	Iterator<TRow>::Type AlignIter1 = begin(row1);
	int i=0;
	while(AlignIter1!=AlignIter){
		if(isGap(AlignIter1))
			::std::cout<<i<<std::endl;
		++i;
		++AlignIter1;

	}
*/

	//typedef String<char> TString;
	//TString str1 = "blablablu";
	//TString str2 ="abab";

	//typedef Align<TString, ArrayGaps> TAlign;
	//TAlign align;
	//resize(rows(align),2);
	//assignSource(row(align,0),str1);
	//assignSource(row(align,1),str2);

	//
	//int score = globalAlignment(align,Score<int>(1,-1,-1,-1),AlignConfig<true,true,false,false>(), Gotoh());

	//::std::cout<<align;

	typedef String<AminoAcid> AminoString;
	AminoString str1 = "PNCFDAKQRTASRPL" ;
	AminoString str2 = "CFDKQKNNRTATRDTA";

	Align<AminoString> ali;
	appendValue(rows(ali),str1);
	appendValue(rows(ali),str2);

	LocalAlignmentFinder<> finder(ali);
	Score<int> scoring(3,-2,-5,-1);

	for(int i= 0;i<3;++i){

		::std::cout<<localAlignment(ali,finder,scoring,0,WatermanEggert())<<std::endl;
		::std::cout<<ali<<std::endl;
	}




	return 0;
}