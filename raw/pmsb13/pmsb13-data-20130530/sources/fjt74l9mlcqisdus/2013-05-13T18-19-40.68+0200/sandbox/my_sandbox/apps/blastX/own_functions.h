// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

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
#include <seqan/score.h>
#include <algorithm>
#include <typeinfo>
#include <seqan/align.h>


using namespace std;
using namespace seqan;

typedef StringSet<String<AminoAcid> > StrSetSA;
typedef Iterator<StringSet<String<AminoAcid> > >::Type TStringSetIterator;
typedef Align<String<AminoAcid>, ArrayGaps> TAlign;

/**
*@brief own_functions.h beinhaltet die Deklaration der Klasse Variable und Match_found.
* Die Klasse Variable speichert alle Benutzereingaben und Match_found speichert Informationen zu den gefundenen matches zwischen den Reads und der Proteindatenbank.
*/
class Variable{
public:
	String<char> fasta_file; /**< Dateipfad zur Proteindatenbank */
	String<char> fastq_file; /**< Dateipfad zu den Reads */
	int seed; /**< Laenge der seeds */
	unsigned size_alp; /**< Auf wieviel Zeichen soll das Alphabet reduziert werden */
	String<unsigned int> numb_alp; /**< Wieviele verschiedene Alphabete soll es geben*/
	int hamming_distance; /**< Groesse der Hamming-Distanz bei nicht exakter Suche */
};


class Match_found{
public:

	String<unsigned int> position_read; /**< Position des Reads im StringSet der uebersetzten Reads */
	String<unsigned int> begin_read; /**< Anfangsposition des matches im Read */
	String<unsigned int> end_read; /**< Endposition des matches im */

	String<unsigned int> position_protein; /**< Position des Proteins im StringSet des uebersetzten Proteins */
	String<unsigned int> begin_protein; /**< Anfangsposition des matches im Protein */
	String<unsigned int> end_protein; /**< Endposition des matches im Protein */

	String<unsigned int> score; /**< Scorevalue von der Verifizierung */


};

template <typename Tchar>
int GET_DATAS(String<char> & file,StringSet<String<Tchar> > & sequence,StringSet<String<char> > & id){ 
	/**
	*@brief GET_DATAS oeffnet die Inputfiles und laedt sie in ihre entsprechenden container
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
//Ueberprueft die Eingabeparameter auf ihre Richtigkeit
int check_values(Variable & comVal,StringSet<String<Dna> > & Reads);
// PARSE_ARGUMENTS.CPP ---------------------------------------------------------------------------



// ALPHABET.CPP ----------------------------------------------------------------------------------
StringSet<String<AminoAcid> > GET_ALPHABET_FORCE();
StringSet<String<AminoAcid> > GET_ALPHABET(unsigned & size_alp, String<unsigned int> & numb_alp);
// ALPHABET.CPP ----------------------------------------------------------------------------------



// FINDER.CPP ------------------------------------------------------------------------------------
int get_read_position(unsigned int & pattern_pos,StringSet<int> & found_reads,int binary);
// SPEICHERT MATCH INFORMATIONEN IN MATCH_FOUND KLASSE
int append_to_match_found(Match_found & seed_found,Pair<unsigned int,unsigned int> & begin_found_prot,Pair<unsigned int,unsigned int> & end_found_prot,Pair<unsigned int,unsigned int> & pattern_found,int & seed,StringSet<int> & found_reads);
// FINDET MATCHES DER LAENGE SEED FUER EIN DATENSATZ VON TRANS_READS UND EINEN DATENSATZ TRANS_PROTEINE UND SPEICHERT
// DIESE IN EINE KLASSE VON MATCH_FOUND
int find_matches(StrSetSA & trans_proteine,StrSetSA & trans_reads,int & seed,int & distance,Match_found & seed_found);
// FUER JEDES ALPHABET WIRD UEBERSETZT, SEED MATCHES GESUCHT, VERIFIZIERT UND ÍN EINER TEXTDATEI AUSGEGEBEN
int FIND_MATCHES_FOR_ALL(StringSet<String<Dna> > & reads,StringSet<String<char> > & readID,StrSetSA & proteine,StringSet<String<char> > & proteinID, int & seed,int & distance,StrSetSA & alphabets);
// FINDER.CPP ------------------------------------------------------------------------------------



// TRANSLATE.CPP ---------------------------------------------------------------------------------
// HASH-FUNKTION GIBT FUER JEDES CODON EIN EINDEUTIGEN INT WIEDER
// DER INT GIBT DIE POSITION IN EINEM ARRAY AN, WELCHER INDIREKT DIE
// AMINOSAEURE GESPEICHERT HAELT
int hash(int x,int y,int z);
// FUNKTION BEINHALTET EIN ARRAY MIT 64 POSITIONEN WELCHE DIE ADRESSE
// EINER AMINOSAEURE BEINHALTEN
int get_Amino_Acid_Pos(int pos);
// BEKOMMT EIN CODON UND GIBT GRUPPENNUMMER DER JEWEILIGEN AMINOSÄURE ZURUECK
int get_translate_from_codon(String<Dna> & actual_triplet,String<AminoAcid> & alphabet,int reduced);
// BEKOMMT EIN READ UND UEBERSETZT ES ZUERST UEBER EINE HASH-FUNKTION IN DIE JEWEILIGE AMINOSAEURE
// DA READING FRAME NICHT BEKANNT IST AUCH IN ALLE 6 MOEGLICHEN READING FRAMES (WIRD NICHT GESPEICHERT)
// SONDERN NUR MIT DER GRUPPEN NUMMER DES JEWEILIGEN ALPHABETES IN TRANS_READS GESPEICHERT 
int translate(StrSetSA & trans_reads,String<Dna> & read,String<AminoAcid> & alphabet,int frame,int reduced);
// UEBERSETZT ALLE READS FUER EIN ALPHABET
int translate_reads(StringSet<String<Dna> > & reads,StrSetSA & trans_reads,String<AminoAcid> & alphabet,int reduced);
// BEKOMMT ALLE PROTEINSEQUENZ UND UEBERSETZT DIESE IN DAS VEREINFACHTE ALPHABET
int translate_database(StrSetSA & trans_proteine,StrSetSA & proteine,String<AminoAcid> & alphabet);
// TRANSLATE.CPP ---------------------------------------------------------------------------------



// VERIFY.CPP ------------------------------------------------------------------------------------
//
int append_to_match_found(Match_found & verify_found, Pair<int,int> & intervall_protein,int position_protein,int &position_read);
//
Pair<int,int> get_position_in_prot(TAlign & align,int & begin,int & end);
//
int known_position(int & begin,int & end,StringSet<Pair<unsigned int, unsigned int> > & pair_position);
//
int verify_seed_match(String<AminoAcid> protein,String<AminoAcid> & read,Match_found & verify_found,int & begin, int & end, int position_protein,int & position_read,int distance);
//
int verify_all(Match_found & seed_found, int & distance,StrSetSA & trans_proteine,StrSetSA & trans_reads,Match_found & verify_found,StringSet<String<Dna> > & reads,StrSetSA & proteine,StrSetSA & alphabets);
// VERIFY.CPP ------------------------------------------------------------------------------------



// OUTPUT.CPP ------------------------------------------------------------------------------------
void write_to_file(Match_found & seed_found, StringSet<String<char> > & proteinID, StringSet<String<char> > & readID,String<unsigned int> numb_alp,int alp);
void write_to_file(Match_found & verify_found, StringSet<String<char> > & proteinID, StringSet<String<char> > & readID,StringSet<String<Dna> > & reads,String<unsigned int> numb_alp,int alp);
// OUTPUT.CPP ------------------------------------------------------------------------------------


#define SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#endif
