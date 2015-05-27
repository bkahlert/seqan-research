#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <iostream>

using namespace seqan;
typedef Index<String<char>, IndexEsa<> > BowTieIndex;
typedef Finder<BowTieIndex> BowTieFinder;
typedef Iterator<BowTieIndex, TopDown< ParentLinks<>>>::Type BowTieIterator;

//global variables
	String<char> text = "ATCAGTATACACAGACAGTATCCAATGCAN";
	String<char> pattern = "ATACA";
	String<char> rText;
	String<char> rPattern;
	String<char> cuttedPattern;
	String<char> rightHalf;
	unsigned textLength = length(text); // == length(rText)
	unsigned patternLength = length(pattern); // == length(rPattern)
	unsigned cutPos = patternLength / 2;
	unsigned oddPattern = patternLength % 2;

//global functions

// returns a string backwards from input
String<char> revert (String<char> word){
	String<char> rWord;
	resize(rWord, length(word));
	unsigned n = length(word);
	for (unsigned i = 0; i< n ; i++){
		rWord[i] = word[n-i-1];
	}
	return rWord;
}

// exactly verifies a string in a text based on a position of a finder (stringmatching)
bool verify (String<char> text , String<char> ver, unsigned pos){
	if ( (pos + length(ver)) > length(text))
			return false;
	bool isVerified = true;
	unsigned i = 0;
	while ((i < length(ver)) && isVerified){
		//std::cout << "pos = " << pos + i <<  " [" << text[pos + i]  << "] =?= [" << ver [i] << "] "<< std::endl;
		isVerified = text[pos + i] == ver[i];
		i++;
	}
	return isVerified;
}

// finds all patterns with no missmatches
void exactSearch (BowTieFinder finder){
	while (find(finder, pattern)) {
	std::cout << pattern << ", found: " << pattern << " at position "<< position(finder) << std::endl;
	}
}

// finds all patterns with one missmatch within the second half of the pattern (but none without missmatch)
void forwardSearch (BowTieFinder finder) {
	for (unsigned i = 0; i < (patternLength- cutPos); i++){

		cuttedPattern= pattern;
		resize (cuttedPattern, cutPos + i);
		//std::cout << cuttedPattern;

		rightHalf = rPattern;
		resize(rightHalf, patternLength - cutPos - i - 1 );
		rightHalf = revert(rightHalf);
		//std::cout << " " << rightHalf << std::endl;

		clear(finder);

		while (find(finder, cuttedPattern)) {

			//std::cout << rightHalf << std::endl;
			//std::cout << position(forwardFinder)<< std::endl;
			if ((position(finder)+ patternLength <= textLength)&&(text[position(finder) + length(cuttedPattern)]!= pattern[i+ cutPos]) && ((length(cuttedPattern)) == (patternLength-1))){
				std::cout << pattern << ", found: " << cuttedPattern << text[position(finder)+ length(cuttedPattern)] << rightHalf << " at position "<< position(finder) << std::endl;
			}else if ((position(finder)+ patternLength <= textLength)&&(text[position(finder) + length(cuttedPattern)] != pattern[i+ cutPos]) && (verify ( text , rightHalf , position(finder) + length(cuttedPattern)+ 1 ))){ 
				std::cout << pattern << ", found: " << cuttedPattern << text[position(finder)+ length(cuttedPattern)] << rightHalf << " at position "<< position(finder) << std::endl;
			}
		}
	}
}

// finds all patterns with one missmatch within the first half of the pattern (but none without missmatch)
void backwardSearch (BowTieFinder finder){
	for (unsigned i = 0; i < cutPos; i++){

		cuttedPattern = rPattern;
		resize (cuttedPattern, patternLength- cutPos + i);
		//std::cout << cuttedPattern;

		rightHalf = pattern;
		resize(rightHalf, cutPos -i - 1);
		rightHalf = revert(rightHalf);
		//std::cout << " " << rightHalf << std::endl;

		clear(finder);

		while (find(finder, cuttedPattern)) {
			if ((position(finder)+ patternLength <= textLength)&&((rText[position(finder) + length(cuttedPattern)]) != rPattern[i + patternLength - cutPos]) && ((length(cuttedPattern)) == (patternLength-1))){
				std::cout << pattern << ", found: " << revert(rightHalf) << rText [position(finder)+ length(cuttedPattern)] << revert(cuttedPattern) << " at position " << textLength - patternLength - position(finder) << std::endl;
			}else if ((position(finder)+ patternLength <= textLength)&&((rText[position(finder) + length(cuttedPattern)]) != rPattern[i + patternLength - cutPos]) && (verify(rText, rightHalf, position(finder) + length(cuttedPattern)+ 1 ))){
				std::cout << pattern << ", found: " << revert(rightHalf) << rText [position(finder)+ length(cuttedPattern)] << revert(cuttedPattern) << " at position " << textLength - patternLength - position(finder) << std::endl;
			}
		}
	}
}

char readPos (String<char> mypattern, unsigned pos){
	return (mypattern[pos]);
}


// char c ist the character that is the substitute at the ith position of the right half of the pattern
void backtrack (BowTieIterator iter, String<char> mypattern, unsigned i, String<char> c, bool isBackward){
	rightHalf = revert(mypattern);
	resize(rightHalf, patternLength - cutPos - i - 1 );
	rightHalf = revert(rightHalf);
	append(c,rightHalf);
	//std::cout << c << std::endl;
	if (goDown(iter,c)){	
			for (unsigned k = 0; k < length(getOccurrences(iter));k++){
				String<char> found = mypattern;
				found[cutPos+i] = c[0];
				if(isBackward){
					std::cout << pattern << ", found: " << revert(found) << " at position "<< textLength - patternLength - getOccurrences(iter)[k] << std::endl;
				}else{
					std::cout << pattern << ", found: " << found << " at position "<< getOccurrences(iter)[k] << std::endl;
				}
			}
	}
}


//main
int main(){
	


	resize(rText, length(text));
	resize(rPattern, patternLength);

	rText = revert(text);
	rPattern = revert(pattern);

	//Forward- and Backward-Index
	BowTieIndex forwardIndex(text);
	BowTieFinder forwardFinder(forwardIndex);
	BowTieIterator forwardIter(forwardIndex);

	BowTieIndex backwardIndex(rText);
	BowTieFinder backwardFinder(backwardIndex);
	BowTieIterator backwardIter(backwardIndex);

	//direct text verification search:
	std::cout << "direct text verification search" << std::endl;
	exactSearch(forwardFinder);
	forwardSearch(forwardFinder);
	backwardSearch(backwardFinder);	

	//backtracking verification search:

	std::cout << std::endl << "backtracking verification search" << std::endl;
	cuttedPattern= pattern;
	resize (cuttedPattern, cutPos);

	//exact Search
	if (goDown(forwardIter,pattern)){
		for (unsigned k = 0; k < length(getOccurrences(forwardIter));k++){;
			std::cout << pattern << ", found: " << pattern << " at position " << getOccurrences(forwardIter)[k] << std::endl;
		}
	}
	
	goRoot(forwardIter);
	goDown(forwardIter, cuttedPattern);
	for (unsigned i = 0; i < (patternLength- cutPos); i++){
		//std::cout << cutPos+i << " = " << readPos(pattern, cutPos+i)<< std::endl;
		if (readPos(pattern,cutPos+i)!= 'A') {
			//std::cout << cutPos+i << " : not 'A'" <<std::endl;
			backtrack(forwardIter, pattern, i, "A", false);
		}
		if (readPos(pattern,cutPos+i)!= 'C') {
			//std::cout << cutPos+i << " : not 'C'" <<std::endl;
			backtrack(forwardIter, pattern, i, "C", false);
		}
		if (readPos(pattern,cutPos+i)!= 'G') {
			//std::cout << cutPos+i << " : not 'G'" <<std::endl;
			backtrack(forwardIter, pattern, i, "G", false);
		}
		if (readPos(pattern,cutPos+i)!= 'T') {
			//std::cout << cutPos+i << " : not 'T'" <<std::endl;
			backtrack(forwardIter, pattern, i, "T", false);
		}
		if (readPos(pattern,cutPos+i)!= 'N') {
			//std::cout << cutPos+i << " : not 'N'" <<std::endl;
			backtrack(forwardIter, pattern, i, "N", false);
		}
		goRoot(forwardIter);
		String<char> newSearch = (String<char>)readPos(pattern,cutPos +i);
		//std::cout << "new Search Letter :" << newSearch << std::endl;
		append(cuttedPattern, newSearch);
		//std::cout << "new Search :" << cuttedPattern << std::endl;
		goDown(forwardIter, cuttedPattern);
	}



	cuttedPattern = rPattern;
	resize (cuttedPattern, patternLength- cutPos);
	//std::cout << "backwardSearch, search for :" << cuttedPattern << std::endl;

	goDown(backwardIter,cuttedPattern);
	for (unsigned i = 0; i < cutPos; i++) {
		//std::cout << patternLength-cutPos+i << " = " << readPos(rPattern,patternLength-cutPos+i)<< std::endl;
		if (readPos(rPattern,patternLength -cutPos+i)!= 'A') {
			//std::cout << patternLength -cutPos+i << " : not 'A'" <<std::endl;
			backtrack(backwardIter, rPattern, i+oddPattern, "A", true);
		}
		if (readPos(rPattern,patternLength -cutPos+i)!= 'C') {
			//std::cout << patternLength -cutPos+i << " : not 'C'" <<std::endl;
			backtrack(backwardIter, rPattern, i+oddPattern, "C", true);
		}
		if (readPos(rPattern,patternLength -cutPos+i)!= 'G') {
			//std::cout << patternLength -cutPos+i << " : not 'G'" <<std::endl;
			backtrack(backwardIter, rPattern, i+oddPattern, "G", true);
		}
		if (readPos(rPattern,patternLength -cutPos+i)!= 'T') {
			//std::cout << patternLength -cutPos+i << " : not 'T'" <<std::endl;
			backtrack(backwardIter, rPattern, i+oddPattern, "T", true);
		}
		if (readPos(rPattern,patternLength -cutPos+i)!= 'N') {
			//std::cout << patternLength -cutPos+i << " : not 'N'" <<std::endl;
			backtrack(backwardIter, rPattern, i+oddPattern, "N", true);
		}
		goRoot(backwardIter);
		String<char> newSearch =(String<char>)readPos(rPattern,patternLength -cutPos +i);
		//std::cout << "new Search Letter :" << newSearch << std::endl;
		append(cuttedPattern, newSearch);
		//std::cout << "new Search :" << cuttedPattern << std::endl;
		goDown(backwardIter, cuttedPattern);

	}
	return 0;
}

