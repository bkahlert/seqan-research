#include <iostream>
#include <seqan/index.h>

using namespace seqan;

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




	typedef Index<DnaString, IndexQGram< OneGappedShape > > TIndex;
	TIndex index("CATGATTACATA");
	stringToShape(indexShape(index),"1101");
	hash(indexShape(index),"AT-T");
	for(unsigned int i =0; i < length(getOccurrences(index, indexShape(index)));++i)
		std::cout<<getOccurrences(index, indexShape(index))[i]<<std::endl;

	

	return 0; 

}