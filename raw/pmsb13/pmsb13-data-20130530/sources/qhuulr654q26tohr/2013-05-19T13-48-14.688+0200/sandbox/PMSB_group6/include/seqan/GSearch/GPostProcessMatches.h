#ifndef GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_
#define GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_

#include <seqan/GStructs/GMatch.h>
#include <seqan/sequence.h>
#include <algorithm>

using namespace seqan;
using namespace std;

bool myVergleich (GMatch<int> & i,GMatch<int> & j) {
	return (i.getScore() >= j.getScore());
}

int sortiereMatches (String<GMatch<int> > & unsortetResults) {
	std::sort (begin(unsortetResults), end(unsortetResults), myVergleich);
	return 0;
}


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
	//quickSort(matches, 0, length(matches)-1);
	//sort (begin(matches), end(matches), myVergleich);
	sortiereMatches(matches);
	return 0;
}

#endif  // GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_