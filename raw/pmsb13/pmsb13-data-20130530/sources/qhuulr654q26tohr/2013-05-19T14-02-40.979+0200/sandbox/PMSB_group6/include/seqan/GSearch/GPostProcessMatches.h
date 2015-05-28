#ifndef GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_
#define GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_

#include <seqan/GStructs/GMatch.h>
#include <seqan/sequence.h>
#include <algorithm>

using namespace seqan;
using namespace std;

template<typename TScore>
void quickSort(String<GMatch<TScore> > & arr, int left, int right) {
	int i = left, j = right;	
	TScore tmp;	
	TScore pivot = arr[(left + right) / 2].getScore();	
	
	/* partition */	
	while (i <= j) {		
		while (arr[i].getScore() < pivot)
			i++;
		while (arr[j].getScore() > pivot)
			j--;
		if (i <= j) {			
			tmp = arr[i].getScore();			
			arr[i].score = arr[j].getScore();			
			arr[j].score = tmp;			
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


template<typename TScore>
int GPostProcessMatches(String<GMatch<TScore> > & matches)
{
	if (length(matches)!=0)
		quickSort(matches, 0, length(matches)-1);
	return 0;
}

#endif  // GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_