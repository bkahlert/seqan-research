// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 18.APRIL.2013

#ifndef SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <cstring>
#include <string>
#include <vector>


using namespace std;
using namespace seqan;

// Klasse Variable speichert alle variablen
// Benutzereingaben
class Variable{
public:
	String<char> fasta_file;
	String<char> fastq_file;
	int seed;
	int size_alp;
	int numb_alp;
};

// OEFFNET INPUTFILES UND LAEDT DIESE IN IHRE ENTSPRECHENDEN CONTAINER
template <typename Tchar>
int getDatas(String<char> file,StringSet<String<Tchar>> & Sequence,StringSet<String<char>> & ID){ 
	SequenceStream seqStream(toCString(file));
	if (!isGood(seqStream)){
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
	readAll(ID, Sequence, seqStream);
	return 0;
}

// PARSE_ARGUMENTS.CPP ---------------------------------------------------------------------------
// INITIALISIERUNG DER KLASSE VARIABLE MIT DEFAULT WERTEN
void dafault_values(Variable & comVal);
// WERTE AUS KOMMANDOZEILE WERDEN GEPARST UND IN comVal GESCHRIEBEN
// WENN DAS PARSEN FEHL SCHLAEGT DANN RETURN WERT 1 
// SONST SCHLIEßt DIE FUNKTION MIT RETURN WERT 0
int PARSE_ARGUMENTS(int argc,char const ** argv,Variable & comVal);
// PARSE_ARGUMENTS.CPP ---------------------------------------------------------------------------



// ALPHABET.CPP ----------------------------------------------------------------------------------
void getAlphabet(int & numb_alp,int & size_alp,StringSet<String<int>> & Alphabete);
// ALPHABET.CPP ----------------------------------------------------------------------------------



// FINDER.CPP ------------------------------------------------------------------------------------
//
void findMatches(StringSet<String<Dna>> & Reads,StringSet<String<char>> & ReadID,StringSet<String<AminoAcid>> & Proteine,StringSet<String<char>> & ProteinID, int seed,StringSet<String<int>> & Alphabete);
// FINDER.CPP ------------------------------------------------------------------------------------


// TRANSLATE.CPP ---------------------------------------------------------------------------------
void translate_reads(StringSet<String<int>> & trans_reads,String<Dna> & read,String<int> & Alphabete,int frame);

int translate(String<Dna> & aktual_triplet,String<int> & Alphabete);
int hash(int x,int y,int z);

// TRANSLATE.CPP ---------------------------------------------------------------------------------

#define SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#endif
