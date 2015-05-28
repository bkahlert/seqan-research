// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 18.APRIL.2013

#ifndef SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <cstring>
#include <string>


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

// PARSE_ARGUMENTS.CPP
void dafault_values(Variable & comVal);
int PARSE_ARGUMENTS(int argc,char const ** argv,Variable & comVal);
// ALPHABET.CPP
void getAlphabet(int & numb_alp,int & size_alp,StringSet<String<int>> & Alphabete);

// OEFFNET INPUTFILES UND LAEDT DIESE IN IHRE ENTSPRECHENDEN
// CONTAINER
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

//
void findMatches(StringSet<String<Dna>> & Reads,StringSet<String<char>> & ReadID,StringSet<String<AminoAcid>> & Proteine,StringSet<String<char>> & ProteinID, int seed,StringSet<String<int>> & Alphabete);
#define SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#endif
