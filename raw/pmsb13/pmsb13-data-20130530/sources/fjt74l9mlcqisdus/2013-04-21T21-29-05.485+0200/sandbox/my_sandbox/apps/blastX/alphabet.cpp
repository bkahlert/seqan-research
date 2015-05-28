// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 18.APRIL.2013

#include "own_functions.h"

// FUNKTION NICHT FERTIG ERSTMAL NUR ZUM WEITERMACHEN SPORALISCH
void getAlphabet(int & numb_alp,int & size_alp,StringSet<String<int>> & Alphabete){
	String<int>alp1; 
	append(alp1,4); // Alanin			0
	append(alp1,3); // Arginin			1
	append(alp1,5); // Asparagin		2
	append(alp1,2); // Asparaginsäure	3
	append(alp1,1); // Cystein			4
	append(alp1,5);	// Glutamin			5
	append(alp1,2); // Glutaminsäure	6
	append(alp1,4); // Glycin			7
	append(alp1,3); // Histidin			8
	append(alp1,4); // Isoleucin		9
	append(alp1,4); // Leucin			10
	append(alp1,3); // Lysin			11
	append(alp1,1); // Methionin		12
	append(alp1,4); // Phenylalanin		13
	append(alp1,4);	// Prolin			14
	append(alp1,5); // Serin			15
	append(alp1,5); // Threonin			16
	append(alp1,5); // Tryptophan		17
	append(alp1,5); // Tyrosin			18
	append(alp1,4); // Valin			19
	append(alp1,6);	// STOP				20
	appendValue(Alphabete,alp1);
}
