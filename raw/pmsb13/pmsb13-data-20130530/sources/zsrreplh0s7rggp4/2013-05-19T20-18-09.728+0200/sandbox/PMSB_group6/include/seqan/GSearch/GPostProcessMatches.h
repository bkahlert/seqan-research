#ifndef GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_
#define GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_

#include <seqan/GStructs/GMatch.h>
#include <seqan/sequence.h>
#include <algorithm>

using namespace seqan;
using namespace std;


// Code Snippet from http://www.algolist.net/Algorithms/Sorting/Quicksort
template<TScore,TSequence>
void quickSort(String<GMatch<TScore,TSequence> > & arr, int left, int right) {
	int i = left, j = right;	
	GMatch<TScore,TSequence> tmp;	
	TScore pivot = arr[(left + right) / 2].getScore();	
	
	/* partition */	
	while (i <= j) {		
		while (arr[i].getScore() > pivot)
			i++;
		while (arr[j].getScore() < pivot)
			j--;
		if (i <= j) {			
			tmp = arr[i];			
			arr[i] = arr[j];			
			arr[j] = tmp;			
			i++;			
			j--;			
		}		
	};	

		/* recursion */

	if (left < j)		
		quickSort(arr, left, j);
		if (i < right)
				quickSort(arr, i, right);
	}


template<TScore,TSequence>
int GPostProcessMatches(String<GMatch<TScore,TSequence> > & matches)
{
	if (length(matches)>=2)
		quickSort(matches, 0, length(matches)-1);
	return 0;
}

#endif  // GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_