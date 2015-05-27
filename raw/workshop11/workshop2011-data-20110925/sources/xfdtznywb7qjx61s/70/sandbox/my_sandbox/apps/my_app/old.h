/*
 *  old.h
 *  seqan
 *
 *  Created by Roland Krause on 13.09.11.
 *  Copyright 2011 MPI for Molecular Genetics. All rights reserved.
 *
 */

void showAllLetterOfMyAlphabet(TAlphabet const &)
{
    typedef typename Size<TAlphabet>::Type TSize;
    TSize alphSize = ValueSize<TAlphabet>::VALUE;
    for (TSize i = 0; i < alphSize; ++i)
        std::cout << i << ',' << TAlphabet(i) << "  ";
    std::cout << std::endl;
}

//template <typename String> class <CharString, DnaString, Peptide>
template <typename TString>
void countOneMers(TString const sequence)
{
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TString>::Type TAlphabet;
	
	typedef String<TSize> TCounterString;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	TCounterString table;
	resize(table, alphSize,0); 
	
	typedef typename Iterator<TString>::Type TIter;
    TIter itEnd = end(sequence);
    for (TIter it = begin(sequence); it != itEnd; goNext(it))
        value(table, ordValue(value(it))) += 1;
	
	
	typedef typename Iterator<TCounterString>::Type TTableIter;	
	TTableIter countIt = begin(table);
    TTableIter countItEnd = end(table);
	for(TSize pos=0; countIt != countItEnd; ++countIt, ++pos)
	{
		if(value(countIt > 0))
		{
			std::cout <<TAlphabet(pos) << ':' << value(countIt) << std::endl;
		}
	}
	
}
