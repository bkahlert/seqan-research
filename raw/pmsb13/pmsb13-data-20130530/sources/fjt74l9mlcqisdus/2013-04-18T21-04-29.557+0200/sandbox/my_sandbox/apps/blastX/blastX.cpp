// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 18.APRIL.2013

#include "own_functions.h"

// MAIN FUNKTION -------------------------------------------------------------------------------
int main(int argc, char const ** argv){
	// Klasse Variable enthält alle Parameter die der
	// Benutzer über Kommandozeile festlegen kann
	Variable comVal;
	// comVal wird mit default Parameter initialisiert
	dafault_values(comVal);
	// Eingabe vom Benutzer wird gelesen und default
	// Parameter wenn Eingabe vorhanden ist ueberschrieben
	if (argc>=2){
		if (PARSE_ARGUMENTS(argc,argv,comVal))
			return 1;
	}
	StringSet<String<int>> x;
	String<int>y;
	append(y,1);
	append(y,2);
	append(y,3);
	appendValue(x,y);

	return 0;
}
// MAIN FUNKTION -------------------------------------------------------------------------------