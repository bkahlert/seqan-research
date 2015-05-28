// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 18.APRIL.2013

#include "own_functions.h"

void getAlphabet(int & numb_alp,int & size_alp,StringSet<String<int>> & Alphabet){
	String<int>alp1; 
	append(alp1,1);
	append(alp1,1);
	append(alp1,2);
	append(alp1,2);
	append(alp1,3);
	append(alp1,3);
	append(alp1,3);
	append(alp1,4);
	append(alp1,4);
	append(alp1,4);
	append(alp1,4);
	append(alp1,4);
	append(alp1,4);
	append(alp1,4);
	append(alp1,5);
	append(alp1,5);
	append(alp1,5);
	append(alp1,5);
	append(alp1,5);
	append(alp1,5);
	appendValue(Alphabet,alp1);
}
