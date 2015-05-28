#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/basic.h>
#include <cstring>
#include <seqan/index_extras.h>
#include <seqan/sequence_journaled.h>
 
using namespace std;
using namespace seqan;

typedef StringSet<String<AminoAcid> > StrSetSA;

int main(){
	StrSetSA trans_proteine;
	StrSetSA trans_reads;
	appendValue (trans_reads,"CRC");
	appendValue (trans_reads,"CQQ");
	appendValue (trans_reads,"RCQ");
	appendValue (trans_reads,"QCR");
	appendValue(trans_proteine,"CRCQQCQQ");
	appendValue(trans_proteine,"QCQQCRCQ");
	Index <StrSetSA, IndexSa<> > index_trans_proteine(trans_proteine);        // hier werden die indexe ¸ber meine ¸bersetzte datenbank und pattern gebaut
    Index <StrSetSA, IndexWotd<> > index_trans_reads(trans_reads);
    
    Finder <Index<StrSetSA, IndexSa<> >,Backtracking<HammingDistance> > finder_obj(index_trans_proteine);
    Pattern <Index<StrSetSA, IndexWotd<> >,Backtracking<HammingDistance> > pattern_obj(index_trans_reads,3);    // laut Jochen muss hier noch die L‰nge des Pattern
    
    while(find(finder_obj,pattern_obj, 1)){     // <- hier kommt es jetzt zu konflikten
    cout <<"["<<beginPosition(finder_obj)<<";"<<endPosition(finder_obj)<<")\t"<<infix(finder_obj)<<"\t"<<position(pattern_obj)<<endl;

    } 
    
    return 0;
}