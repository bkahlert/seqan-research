#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <iostream>

using namespace seqan;

///*
//template <typename T>
//struct Iterator;
//
//template <typename TElement>
//struct Iterator<String<TElement> >
//{
//    typedef Iter<String<TElement>, StringIterator> Type;
//}
//
//Iterator<String<TElement> >::Type;
//
//class A {
//  class B {
//    typedef int X;
//  };
//};
//
//A::B::X y;
//*/
//
//String<AminoAcid> sequenz = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
//Iterator<String<AminoAcid> >::Type ende = end(sequenz);
//
//
//int main() {
//	
//	for ( Iterator<String<AminoAcid> >::Type it = begin(sequenz) ; it != ende ; goNext(it))
//	{
//		if (value(it) == 'R')
//			value(it) = 'A';
//	}
//	
//	std::cout << sequenz << std::endl << std::endl;
//
//	String<Size<String<AminoAcid> >::Type> counter;
//	Size<String<AminoAcid> >::Type alphSize = ValueSize<AminoAcid>::VALUE;
//
//	std::cout << ValueSize<AminoAcid>::VALUE << " " << alphSize << std::endl;
//
//	resize(counter, alphSize, 0);
//
//	std::cout << ValueSize<AminoAcid>::VALUE << " " << alphSize << std::endl;
//
//	return 0;
//
//}

 int main() {
	String<char> Text = "GTATACACAGATAGTß";
	String<char> Pattern = "ATACA";
 
	Index<String<char>, IndexEsa<> > index;
	std::cout << "bla2" << std::endl;
 
 	return 0;
 }
