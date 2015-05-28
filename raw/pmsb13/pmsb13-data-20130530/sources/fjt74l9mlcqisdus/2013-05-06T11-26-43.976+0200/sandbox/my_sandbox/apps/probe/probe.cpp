#include <iostream>
#include <seqan/index.h>
#include <iostream>
#include <seqan/index_extras.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

typedef StringSet<String<AminoAcid>> StrSetSA; 

using namespace std;

using namespace seqan;

int main(){
	StrSetSA trans_proteine;
	StrSetSA trans_reads;
	appendValue (trans_reads,"CRC");
	appendValue (trans_reads,"CQQ");
	appendValue(trans_proteine,"CRCQQCQQ");
	appendValue(trans_proteine,"QCQQCRCQ");
	
	Index <StrSetSA> index_trans_proteine(trans_proteine);        // hier werden die indexe �ber meine �bersetzte datenbank und pattern gebaut
    Index <StrSetSA> index_trans_reads(trans_reads);
    
    Finder <Index<StrSetSA>,Backtracking<HammingDistance> > finder_obj(index_trans_proteine);
    Pattern <Index<StrSetSA>,Backtracking<HammingDistance> > pattern_obj(index_trans_reads,length(trans_reads));    // laut Jochen muss hier noch die L�nge des Pattern
                                                                                                                                                                                     �bergeben werden
        
    while(find(finder_obj,pattern_obj,0)){     // <- hier kommt es jetzt zu konflikten
    cout <<"["<<beginPosition(finder_obj)<<";"<<endPosition(finder_obj)<<")\t"<<infix(finder_obj)<<endl;


}