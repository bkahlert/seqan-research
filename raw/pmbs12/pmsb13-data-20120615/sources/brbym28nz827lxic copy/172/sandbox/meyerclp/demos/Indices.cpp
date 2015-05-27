#include <iostream>
#include <seqan/index.h>
#include <seqan/align.h>

using namespace seqan;



template <typename TStringSet, typename TIndexSpec>
void qgramCounting(TStringSet &set, TIndexSpec)
{
    typedef Index<TStringSet, TIndexSpec> TIndex;
    typedef typename Fibre<TIndex, QGramCounts>::Type TCounts;
    typedef typename Fibre<TIndex, QGramCountsDir>::Type TCountsDir;
    typedef typename Value<TCountsDir>::Type TDirValue;
    typedef typename Iterator<TCounts, Standard>::Type TIterCounts;
    typedef typename Iterator<TCountsDir, Standard>::Type TIterCountsDir;
    TIndex index(set);
    indexRequire(index, QGramCounts());
    // initialize distance matrix
    int seqNum = countSequences(index);
    Matrix<int, 2> distMat;
    setLength(distMat, 0, seqNum);
    setLength(distMat, 1, seqNum);
    resize(distMat, 0);
    std::cout << std::endl << "Length of the CountsDir fibre: " << length(indexCountsDir(index)) << std::endl;
    TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
    TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
    TIterCounts itCountsBegin = begin(indexCounts(index), Standard());

	 // for each bucket count common q-grams for each sequence pair
    TDirValue bucketBegin = *itCountsDir;
    for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir)
    {
        TDirValue bucketEnd = *itCountsDir;
        // q-gram must occur in at least 2 different sequences
        if (bucketBegin != bucketEnd)
        {
            TIterCounts itA = itCountsBegin + bucketBegin;
            TIterCounts itEnd = itCountsBegin + bucketEnd;
            for(; itA != itEnd; ++itA)
                for(TIterCounts itB = itA; itB != itEnd; ++itB)
                    distMat((*itA).i1, (*itB).i1) += _min((*itA).i2, (*itB).i2);
        }
        bucketBegin = bucketEnd;
    }
    std::cout << std::endl << "Common 5-mers for Seq_i, Seq_j" << std::endl;
    std::cout << distMat;
}






int main (){

	//String<char> myString = "abracadabra";
 //   typedef Index< String<char> > TMyIndex;
 //   TMyIndex myIndex(myString);

	//Iterator< TMyIndex, TopDown< ParentLinks<Preorder> > >::Type myIterator(myIndex);
 //   while (!atEnd(myIterator))
 //   {
 //       std::cout << representative(myIterator) << std::endl;
 //       ++myIterator;
 //   }

	//typedef String<char> TString;
	//StringSet<TString> seq;
	//appendValue(seq,"tobeornottobe");
	//appendValue(seq,"thebeeonthecomb");
	//appendValue(seq,"beingjohnmalkovich");

	//typedef Index< StringSet<TString> > TMyIndex;
	//TMyIndex myIndex(seq);

	//Iterator< TMyIndex, BottomUp<>>::Type myIterator(myIndex);
	//while (!atEnd(myIterator))
 //   {
 //       std::cout << representative(myIterator) << std::endl;
 //       ++myIterator;
 //   }
	
	//typedef String<char> TString;
	//StringSet<TString> seq;
	//appendValue(seq,"CDFGHC");
	//appendValue(seq,"CDEFGAHC");

	//typedef Index< StringSet<TString> > TMyIndex;
	//TMyIndex myIndex(seq);

	//Iterator< TMyIndex, Mums >::Type myIterator(myIndex);
	//	while (!atEnd(myIterator))
 //   {
 //       std::cout << representative(myIterator) << std::endl;
 //       ++myIterator;
 //   }




	//typedef String<char> TString;
	//TString seq = "tobeornottobe";

	//typedef Index<TString> TIndex;
	//TIndex myIndex(seq);
	//Iterator<TIndex, TopDown<ParentLinks<> > >::Type it(myIndex);
	//bool a= true;
	//while(a){
	//	
	//	
	//	if(goDown(it)){
	//		::std::cout<<representative(it)<<std::endl;
	//		
	//	}
	//	else{
	//		if(goRight(it))
	//			::std::cout<<representative(it)<<std::endl;
	//		else{
	//			while(goUp(it)){
	//				if(goRight(it)){
	//					::std::cout<<representative(it)<<std::endl;
	//					break;
	//				}
	//				a=false;

	//			}
	//		}
	//	}
	//}




	//typedef Index<DnaString, IndexQGram< OneGappedShape > > TIndex;//da gap erst noch stringToShape
	//TIndex index("CATGATTACATA");
	//stringToShape(indexShape(index),"1101");
	//hash(indexShape(index),"AT-A");
	//for(unsigned int i =0; i < length(getOccurrences(index, indexShape(index)));++i)
	//	std::cout<<getOccurrences(index, indexShape(index))[i]<<std::endl;
	
	
	//  for the sake of reproducibility
    Rng<MersenneTwister> rng;
    // create StringSet of 3 random sequences
    StringSet<DnaString> stringSet;
    reserve(stringSet, 3);
    for (int seqNo = 0; seqNo < 3; ++seqNo)
    {
        DnaString tmp;
        int len = pickRandomNumber(rng) % 100 + 10;
        for (int i = 0; i < len; ++i)
            appendValue(tmp, Dna(pickRandomNumber(rng) % 4));
        appendValue(stringSet, tmp);
        std::cout << ">Seq" << seqNo << std::endl << tmp << std::endl;
    }
    qgramCounting(stringSet, IndexQGram<UngappedShape<5> >());
    qgramCounting(stringSet, IndexQGram<UngappedShape<5>, OpenAddressing>());
	



	return 0; 

}