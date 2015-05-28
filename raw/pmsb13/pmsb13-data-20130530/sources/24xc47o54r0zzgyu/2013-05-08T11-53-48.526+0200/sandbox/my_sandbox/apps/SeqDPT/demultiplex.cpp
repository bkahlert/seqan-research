///*
///*! Program part for barcode demultiplexing */
//
//#include <seqan/find.h>
//#include <seqan/sequence.h>
//#include <seqan/index.h>
//#include <seqan/seq_io.h>
//
//using namespace seqan;
//
//
//typedef Dna5String TAlphabet;
////Index<TAlphabet> index;
////Finder<Index<CharString> > finder;
////StringSet<TAlphabet> needles;
//
//
///*!
//* Function for creating the haystacks of barcodes.
//*/
////template <typename TBcs, typename TIndices, typename TFinders> //must get an empty vector for indices and one for finder
////void makeExactIndices(const TBcs& barcodes, TIndices& indices, TFinders& finders)
////{
////	typedef Iterator<Tbcs>::Type TIterator;
////	for (TIterator it = begin(bcs); it != end(bcs); goNext(it))
////	{
////		Index<TAlphabet, IndexEsa<> > esaIndex(it); 
////		Finder<esaIndex::Type, IndexEsa<> > > esaFinder(esaIndex);
////		append(indices, esaIndex);
////		append(finders, esaFinder);
////	}
////}
//template <typename TBcs, typename TIndices, typename TFinders> 
//void makeExactIndices(const TBcs& barcodes, TIndices& indices, TFinders finders)
//{
//	Index<TAlphabet,IndexEsa<> > esaIndex(barcodes);
//	Finder<esaIndex::Type, IndexEsa<> > esaFinder(esaIndex);
//}
//
///*!
//* Function for searching 1 piece of sequence in the barcode indices
//*/
//template <typename TSeq, typename TFinders>							//must get the vector/array of finders
//int findExactIndex(const TSeq& seq, TFinders& finders)
//{
//	typedef Iterator<TFinders, Rooted>::Type TIterator;
//	for(TIterator it = begin(finders); it != end(finders); goNext(it))
//	{
//		if(find(it, seq)) 
//		{
//			return position(it);									//returns index of barcode -> ONLY THE FIRST HIT!
//		}
//		clear(it);
//	}
//	return -1;														//return -1 if no hit occured
//}
//
///*!
//* Function for performing the search on all given sequence pieces and indices
//*/
//template <typename TSeqs, typename TFinders>
//std::vector<std::pair<unsigned, int> > findAllExactIndex(const TSeqs& seqs, TFinders& finders) 
//{
//	std::vector<pair<unsigned, int> > matches (length(seqs));		//initialises the result vector
//	typedef Iterator<TSeqs, Rooted>::Type TIterator;
//	for(TIterator it = begin(seqs); it != end(seqs); goNext(it))
//	{
//		append(matches, std::pair(position(it), findExactIndex(it, Finders)));
//	}
//	return matches;
//}

