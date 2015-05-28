// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 18.APRIL.2013

#include "own_functions.h"

// FUNKTION NICHT FERTIG ERSTMAL NUR ZUM WEITERMACHEN SPORALISCH
void getAlphabet(int & numb_alp,int & size_alp,StringSet<String<int>> & Alphabete){
	String<int>alp1; 
	append(alp1,1); // Methionin
	append(alp1,1); // Cystein
	append(alp1,2); // Glutaminsäure
	append(alp1,2); // Asparaginsäure
	append(alp1,3); // Lysin 
	append(alp1,3);	// Arginin
	append(alp1,3); // Histidin
	append(alp1,4); // Alanin
	append(alp1,4); // Prolin
	append(alp1,4); // Valin
	append(alp1,4); // Leucin
	append(alp1,4); // Isoleucin
	append(alp1,4); // Glycin
	append(alp1,4); // Phenylalanin
	append(alp1,5);	// Serin
	append(alp1,5); // Asparagin
	append(alp1,5); // Glutamin
	append(alp1,5); // Threonin
	append(alp1,5); // Tyrosin
	append(alp1,5); // Typtophan
	append(alp1,6);	// STOP
	appendValue(Alphabete,alp1);
}
