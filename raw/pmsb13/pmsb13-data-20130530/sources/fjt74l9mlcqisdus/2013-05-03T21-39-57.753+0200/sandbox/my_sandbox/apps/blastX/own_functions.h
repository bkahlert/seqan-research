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
#include <seqan/basic.h>
#include <cstring>
#include <seqan/index_extras.h>

using namespace std;
using namespace seqan;

typedef StringSet<String<AminoAcid> > StrSetSA;

/**
*@brief In dieser Datei befindet sich die Deklaration der Klasse Variable.
* Die Klasse speichert alle Benutzereingaben.
*/
class Variable{
public:
	String<char> fasta_file; /**< Dateipfad zur Proteindatenbank */
	String<char> fastq_file; /**< Dateipfad zu den Reads */
	int seed; /**< Laenge der seeds */
	int size_alp; /**< Auf wieviel Zeichen soll das Alphabet reduziert werden */
	int numb_alp; /**< Wieviele verschiedene Alphabete soll es geben*/
	int hamming_distance; /**< Groesse der Hamming-Distanz bei nicht exakter Suche */
};


class Match_found{
public:

	StringSet<unsigned int> position_read; /**< Position Read in StringSet trans_reads */
	StringSet<unsigned int> begin_read; /**< begin in read */
	StringSet<unsigned int> end_read; /**< end in read */

	StringSet<unsigned int> position_protein; /**< begin in protein */
	StringSet<unsigned int> begin_protein; /**< begin in protein */
	StringSet<unsigned int> end_protein; /**< end in protein */
};

template <typename Tchar>
int GET_DATAS(String<char> file,StringSet<String<Tchar> > & sequence,StringSet<String<char> > & id){ 
	/**
	*@brief Oeffnet die Inputfiles und laedt sie in ihre entsprechenden container
	*@param file Name der Inputfile
	*@param sequence Sequenzen in den Files
	*@param id Ids der verschiedenen Sequenzen
	*@return Gibt 0 zurueck wenn die Daten erfolgreich gespeichert wurden und 1 falls es fehlschlaegt
	*/
	SequenceStream seqStream(toCString(file));
	if (!isGood(seqStream)){
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
	readAll(id, sequence, seqStream);
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
StrSetSA GET_ALPHABET(int & numb_alp,int & size_alp);
// ALPHABET.CPP ----------------------------------------------------------------------------------



// FINDER.CPP ------------------------------------------------------------------------------------
//
int find_matches(StrSetSA & trans_proteine,StrSetSA & trans_reads,int seed,Match_found seed_found);
// 
int FIND_MATCHES_FOR_ALL(StringSet<String<Dna5> > & Reads,StringSet<String<char> > & ReadID,StrSetSA & Proteine,StringSet<String<char> > & ProteinID, int seed,StrSetSA & Alphabete);
// FINDER.CPP ------------------------------------------------------------------------------------



// TRANSLATE.CPP ---------------------------------------------------------------------------------
// BEKOMMT EIN READ UND UEBERSETZT ES ZUERST UEBER EINE HASH-FUNKTION IN DIE JEWEILIGE AMINOSAEURE
// DA READING FRAME NICHT BEKANNT IST AUCH IN ALLE 6 MOEGLICHEN READING FRAMES (WIRD NICHT GESPEICHERT)
// SONDERN NUR MIT DER GRUPPEN NUMMER DES JEWEILIGEN ALPHABETES IN TRANS_READS GESPEICHERT 
int translate_reads(StrSetSA & trans_reads,String<Dna5> & read,String<AminoAcid> & Alphabete,int frame);
// BEKOMMT EIN CODON UND GIBT GRUPPENNUMMER DER JEWEILIGEN AMINOSÄURE ZURUECK
int get_translate(String<Dna5> & aktual_triplet,String<AminoAcid> & Alphabete);
// HASH-FUNKTION GIBT FUER JEDES CODON EIN EINDEUTIGEN INT WIEDER
// DER INT GIBT DIE POSITION IN EINEM ARRAY AN, WELCHER INDIREKT DIE
// AMINOSAEURE GESPEICHERT HAELT
int hash(int x,int y,int z);
// FUNKTION BEINHALTET EIN ARRAY MIT 64 POSITIONEN WELCHE DIE ADRESSE
// EINER AMINOSAEURE BEINHALTEN
int get_Amino_Acid_Pos(int pos);
// BEKOMMT ALLE PROTEINSEQUENZ UND UEBERSETZT DIESE IN DAS VEREINFACHTE ALPHABET
int translate_datenbank(StrSetSA & trans_proteine,StrSetSA & protein,String<AminoAcid> & Alphabete);
// UEBERSETZT ALLE READS FUER EIN ALPHABET
int get_translate_for_all(StringSet<String<Dna5>> & Reads,StrSetSA & trans_reads,String<AminoAcid> & Alphabete);
// TRANSLATE.CPP ---------------------------------------------------------------------------------

#define SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#endif
