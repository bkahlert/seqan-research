// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 18.APRIL.2013

#ifndef SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#include <iostream>
#include <seqan/arg_parse.h>

using namespace std;
using namespace seqan;

class Variable{
public:
	String<char> fasta_file;
	String<char> fastq_file;
	int seed;
	int size_alp;
	int numb_alp;
};

void dafault_values(Variable & comVal);
int PARSE_ARGUMENTS(int argc,char const ** argv,Variable & comVal);
void getAlphabet(int & numb_alp,int & size_alp,StringSet<String<int>> & Alphabete);
template <typename Tchar>
void getData(String<char> file,StringSet<String<Tchar>> & Sequence,StringSet<String<char>> & ID);


#define SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#endif
