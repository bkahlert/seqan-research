// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 21.APRIL.2013

#ifndef SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <cstring>
#include <string>
#include <vector>


#include <seqan/graph_types.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>

using namespace std;
using namespace seqan;

typedef StringSet<String<int>> StrSetSInt;

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
int GET_DATAS(String<char> file,StringSet<String<Tchar>> & Sequence,StringSet<String<char>> & ID){ 
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
void DEFAULT_VALUES(Variable & comVal);
// WERTE AUS KOMMANDOZEILE WERDEN GEPARST UND IN comVal GESCHRIEBEN
// WENN DAS PARSEN FEHL SCHLAEGT DANN RETURN WERT 1 
// SONST SCHLIEßt DIE FUNKTION MIT RETURN WERT 0
int PARSE_ARGUMENTS(int argc,char const ** argv,Variable & comVal);
// PARSE_ARGUMENTS.CPP ---------------------------------------------------------------------------



// ALPHABET.CPP ----------------------------------------------------------------------------------
void getAlphabet(int & numb_alp,int & size_alp,StrSetSInt & Alphabete);
// ALPHABET.CPP ----------------------------------------------------------------------------------



// FINDER.CPP ------------------------------------------------------------------------------------
//
void FIND_MATCHES(StringSet<String<Dna>> & Reads,StringSet<String<char>> & ReadID,StringSet<String<AminoAcid>> & Proteine,StringSet<String<char>> & ProteinID, int seed,StrSetSInt & Alphabete);
// FINDER.CPP ------------------------------------------------------------------------------------



// TRANSLATE.CPP ---------------------------------------------------------------------------------
// BEKOMMT EIN READ UND UEBERSETZT ES ZUERST UEBER EINE HASH-FUNKTION IN DIE JEWEILIGE AMINOSAEURE
// DA READING FRAME NICHT BEKANNT IST AUCH IN ALLE 6 MOEGLICHEN READING FRAMES (WIRD NICHT GESPEICHERT)
// SONDERN NUR MIT DER GRUPPEN NUMMER DES JEWEILIGEN ALPHABETES IN TRANS_READS GESPEICHERT 
void translate_reads(StrSetSInt & trans_reads,String<Dna> & read,String<int> & Alphabete,int frame);
// BEKOMMT EIN CODON UND GIBT GRUPPENNUMMER DER JEWEILIGEN AMINOSÄURE ZURUECK
int get_translate(String<Dna> & aktual_triplet,String<int> & Alphabete);
// HASH-FUNKTION GIBT FUER JEDES CODON EIN EINDEUTIGEN INT WIEDER
// DER INT GIBT DIE POSITION IN EINEM ARRAY AN, WELCHER INDIREKT DIE
// AMINOSAEURE GESPEICHERT HAELT
int hash(int x,int y,int z);
// FUNKTION BEINHALTET EIN ARRAY MIT 64 POSITIONEN WELCHE DIE ADRESSE
// EINER AMINOSAEURE BEINHALTEN
int get_AminoAcidPos(int pos);
// BEKOMMT EINE PROTEINSEQUENZ UND UEBERSETZT DIESE IN DAS VEREINFACHTE ALPHABET
void translate_datenbank(StringSet<String<int>> & trans_proteine,String<AminoAcid> & protein,String<int> & Alphabete);
// TRANSLATE.CPP ---------------------------------------------------------------------------------

#define SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#endif
