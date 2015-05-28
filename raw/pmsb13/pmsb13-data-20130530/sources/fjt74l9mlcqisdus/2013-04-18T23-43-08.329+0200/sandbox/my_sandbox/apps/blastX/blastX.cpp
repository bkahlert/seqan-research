// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 18.APRIL.2013

#include "own_functions.h"

template <typename Tchar>
int getDatas(String<char> file,StringSet<String<Tchar>> & Sequence,StringSet<String<char>> & ID){ 
	SequenceStream seqStream(toCString(file));
	if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
	//readRecord(ID, Sequence, seqStream);
	return 0;
}



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
	// Alphabet enthält nach getAlphabet numb_alp mal 
	// String<int> der laenge zwanzig wobei jede Position
	// eine AS repraesentiert und der Eintrag der  
	// jeweiligen zugehoerigen Gruppe (reduziertes Alphabet)
	// in diesem Alphabet
	StringSet<String<int>> Alphabete;
	getAlphabet(comVal.numb_alp,comVal.size_alp,Alphabete);
	
	StringSet<String<Dna>> Reads;
	StringSet<String<char>> ReadID;
	getDatas(comVal.fastq_file,Reads,ReadID);

	StringSet<String<AminoAcid>> Proteine;
	StringSet<String<char>> ProteinID;
	//getData(comVal.fasta_file,Proteine,ProteinID);
	return 0;
}
// MAIN FUNKTION -------------------------------------------------------------------------------