/*! Program part for barcode demultiplexing */

#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace seqan;


typedef Dna5String TAlphabet;
//Index<TAlphabet> index;
//Finder<Index<CharString> > finder;
//StringSet<TAlphabet> needles;


template <typename TIterator, typename TPosition>
struct occurence
{
	Iterator<StringSet<TAlphabet> > seq;
	std::pair<TIterator, TPosition> match;
};

/*!
* Function for creating the haystacks of barcodes.
*/
template <typename TBcs, typename TIndices, typename TFinders> //must get an empty vector/array for indices and one for finder
void makeExactIndices(const TBcs& barcodes, TIndices& indices, TFinders& finders)
{
	typedef Iterator<Tbcs>::Type TIterator;
	for (TIterator it = begin(bcs); it != end(bcs); goNext(it))
	{
		Index<TAlphabet, IndexEsa<> > esaIndex(it); 
		Finder<esaIndex::Type, IndexEsa<> > > esaFinder(esaIndex);
		append(indices, esaIndex);
		append(finders, esaFinder);
	}
}

/*!
* Function for searching 1 piece of sequence in the barcode indices
*/
template <typename TSeq, typename TFinders>							//must get the vector/array of finders
std::pair<Iterator<TFinders>, Position<Iterator<TFinders> > > findExactIndex(const TSeq& seq, TFinders& finders)
{
	typedef Iterator<TFinders>::Type TIterator;
	for(TIterator it = begin(finders); it != end(finders); goNext(it))
	{
		if(find(it, seq)) return (std::pair(it, position(it)));		//returns pointer at barcode and position (pair)  -> ONLY THE FIRST HIT!
	}
	return (NULL, NULL);											//Null-pointer if no matches occured
}
/*!
* Function for performing the search on all given sequence pieces and indices
*/
template <typename TSeqs, typename TFinders>
std::vector<occurence<Iterator<TSeqs>, std::pair<Iterator<TFinders>, Position<Iterator<TFinders> > > > > findAllExactIndex(const TSeqs& seqs, TFinders& finders)
{
	std::vector<occurence<Iterator<TSeqs>, std::pair<Iterator<TFinders>, Position<Iterator<TFinders> > > > > matches (length(seqs)); //initialises the result vector
	typedef Iterator<TSeqs>::Type TIterator;
	for(TIterator it = begin(seqs); it != end(seqs); goNext(it))
	{
		std::pair<Iterator<TFinders>, Position<Iterator<TFinders> res = findExactIndex(it, Finders);
		if res != std::pair(NULL, NULL)								//Checkt if the sequence piece could be matched with a barcode
		{
			append(matches, occurence(it, findExactIndex(it, Finders)));
		}
	}
	return matches;
}
