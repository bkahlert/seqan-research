#include <iostream>
#include <seqan/index.h>
#include <iostream>

typedef StringSet<String<AminoAcid>>StrSetSA; 

int main(){
	StrSetSA trans_proteine;
	StrSetSA trans_reads;
	appendValue (trans_reads,"CRC");
	appendValue (trans_reads,"CQQ");
	appendValue(trans_proteine,"CRCQQCQQ");
	appendValue(trans_proteine,"QCQQCRCQ");
	
	Index <StrSetSA> index_trans_proteine(trans_proteine);        // hier werden die indexe über meine übersetzte datenbank und pattern gebaut
    Index <StrSetSA> index_trans_reads(trans_reads);
    
    Finder <Index<StrSetSA>,Backtracking<HammingDistance> > finder_obj(index_trans_proteine);
    Pattern <Index<StrSetSA>,Backtracking<HammingDistance> > pattern_obj(index_trans_reads,length(trans_reads));    // laut Jochen muss hier noch die Länge des Pattern
                                                                                                                                                                                     übergeben werden
        
    while(find(finder_obj,pattern_obj,0)){     // <- hier kommt es jetzt zu konflikten
    cout <<"["<<beginPosition(finder_obj)<<";"<<endPosition(finder_obj)<<")\t"<<infix(finder_obj)<<endl;


}