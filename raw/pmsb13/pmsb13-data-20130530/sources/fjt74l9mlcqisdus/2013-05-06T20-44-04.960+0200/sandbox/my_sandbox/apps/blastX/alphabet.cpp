// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

/**
*@brief alphabet.cpp
* Beinaltet Funktionen, die verschiedene Alphabete erstellen.
*/

#include "own_functions.h"



StringSet<String<AminoAcid> > GET_ALPHABET_FORCE()
{
	/**@brief GET_ALPHABET_FORCE gibt ein reduziertes Aminosaeurealphabet zurueck.
	* Das Alphabet enthaelt sechs Zeichen, wobei die Aminosaeuren
	* nach ihren Kraeften gruppiert sind, die sie in der Sekundaer -bzw.  
	* Tertiaerstruktur ausbilden (z.B: Wasserstoffbruecken, Disulfidbruecken, ...).
	*@return Gibt einen String mit dem reduzierten Alphabet zurueck
	*/
	StringSet<String<AminoAcid> > alphabets;
	String<AminoAcid> alphabet; 
	append(alphabet,4); // Alanin			0
	append(alphabet,3); // Arginin			1
	append(alphabet,5); // Asparagin		2
	append(alphabet,2); // Asparaginsäure	3
	append(alphabet,1); // Cystein			4
	append(alphabet,5);	// Glutamin			5
	append(alphabet,2); // Glutaminsäure	6
	append(alphabet,4); // Glycin			7
	append(alphabet,3); // Histidin			8
	append(alphabet,4); // Isoleucin		9
	append(alphabet,4); // Leucin			10
	append(alphabet,3); // Lysin			11
	append(alphabet,1); // Methionin		12
	append(alphabet,4); // Phenylalanin		13
	append(alphabet,4);	// Prolin			14
	append(alphabet,5); // Serin			15
	append(alphabet,5); // Threonin			16
	append(alphabet,5); // Tryptophan		17
	append(alphabet,5); // Tyrosin			18
	append(alphabet,4); // Valin			19
	append(alphabet,6);	// STOP				20
	appendValue(alphabets,alphabet);
	return (alphabets);
}

String<AminoAcid> get_alphabet(int & score_matrix){
	String<AminoAcid> alphabet;
	resize(alphabet, 21);
	
	//Referenzstring von Aminosaeuren dient zum befuellen der scoreMatrix
	String<AminoAcid> aminoAcid;
	for(unsigned i = 0; i < 20; ++i)	append(aminoAcid, i);
		
	//Anzahl an Clusterdurchlaeufen
	unsigned run = 20 - alpGroup;						
	//pos1 und pos2 speichern die Positionen des Maximums der scoreMatrix	
	int pos1; 
	int pos2;
	//Nummer der Gruppe	
	int count = 0;	
	int scoreMatrix [20][20];
		
	switch(score_matrix){
	case 0:
		Score<int, ScoreMatrix<AminoAcid, Blosum30_> > scoreOfMatrix;
		break;
	case 1:
		Score<int, ScoreMatrix<AminoAcid, Blosum60_> > scoreOfMatrix;
		break;
	case 2:
		Score<int, ScoreMatrix<AminoAcid, Blosum80_> > scoreOfMatrix;
		break;
	case 3:
		Score<int, ScoreMatrix<AminoAcid, Pam120_> > scoreOfMatrix;
		break;
	case 4:
		Score<int, ScoreMatrix<AminoAcid, Pam200_> > scoreOfMatrix;
		break;
	case 5:
		Score<int, ScoreMatrix<AminoAcid, Pam250_> > scoreOfMatrix;
		break;
	case 6:
		Score<int, ScoreMatrix<AminoAcid, Pam40_> > scoreOfMatrix;
		break;
	case 7:
		Score<int, ScoreMatrix<AminoAcid, Vtml200_> > scoreOfMatrix;
		break;
	case 8:
		Score<int, ScoreMatrix<AminoAcid, Blosum45> > scoreOfMatrix;
		break;
	}
	
	//scoreMatrix wird mit Werten befuellt 	
	for(unsigned i = 0; i < length(aminoAcid); ++i) 
	{
		for(unsigned j = 0; j < length(aminoAcid); ++j)
		{
			scoreMatrix[i][j] = score(scoreOfMatrix, aminoAcid[i], aminoAcid[j]);
			if(i == j)
				scoreMatrix[i][j] = -50;	
		}
	}
	for(unsigned r = 0; r < run; ++r)
	{
		//Speichert das Maximum der score Matrix
		int maxWert = -50;	
		for(unsigned k = 0; k < 20; ++k)
		{		
			for(unsigned l = 0; l < 19; ++l)
			{
				if(max(scoreMatrix[k][l], scoreMatrix[k][l+1]) > maxWert)			
				{	
					maxWert = max(scoreMatrix[k][l], scoreMatrix[k][l+1]);
				}			
			}
		}	 
		//Sucht die Position des Maximums
		for(unsigned i = 0; i < 20; ++i) 
		{
			for(unsigned j = 0; j < 20; ++j)
			{
				if(scoreMatrix[i][j] == maxWert)
				{	
					pos1 = i;
					pos2 = j;
				}
			}
		}
		//Setzt score auf -50, damit er beim naechsten Durchgang nicht wieder gefunden wird			
		scoreMatrix[pos1][pos2] = -50; 
		scoreMatrix[pos2][pos1] = -50;
		//Wenn beide AS noch keiner Gruppe zugeordnet sind		
		if(!alphabet[pos1] && !alphabet[pos2])
		{					
			++count;
			alphabet[pos1] = count;
			alphabet[pos2] = count;
		}
		//Wenn einer von den AS bereits einer Gruppe angehoert uebernimmt die andere diese...
		if(alphabet[pos1] && !alphabet[pos2])
		{
			alphabet[pos2] = alphabet[pos1];
		}	
		//...und andersherum
		if(!alphabet[pos1] && alphabet[pos2])
		{
			alphabet[pos1] = alphabet[pos2];
		}
	}
	//Teilt die restlichen AS gruppen zu
	for(unsigned i = 0; i < length(alphabet); ++i)
	{
		if(!alphabet[i])		
		{	
			++count;	
			alphabet[i] = count;
		}
	}
	return(alphabet);
}

StringSet<String<AminoAcid> > GET_ALL_ALPHABETS(int & alpGroup, int & varianten)
{
	/**@brief GET_ALPHABET gibt ein reduziertes Aminosaeurealphabet zurueck.
	*Alphabet entseht durch Clustern mit zur Hilfenahme von Scorematrizen (z.B. Blosum62, PAM120, etc..)
	*@param matrixNumb gibt an wieviel verschiedene Scorematrizen benutzt werden 
	*@param alpGroup gibt an auf wieviel Cluster das Alphabet reduziert werden soll plus Stoppcodon
	*@return Gibt einen String mit reduzierten Alphabet zurueck
	*/
		
	StringSet<String<AminoAcid> > alphabets;

	//Wird mit den Gruppennummern der jeweiligen AS befuellt, Position im String entspricht AS	
	for (unsigned var=0;var<varianten;++var){
		String<AminoAcid> alphabet = get_alphabet(seqanScorematrices[var]);
		appendValue(alphabets,alphabet);
	}
		
	return (alphabets);
}
