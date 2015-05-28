#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc,char * argv[])
{
   if(argc<2){
   std::cerr<<"Ungueltige Anzahl an Parametern"<<std::endl;
   return -1;
   }
    char * file=argv[1];
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::SequenceStream seqStream(file);
    if(seqStream.good()){
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << std::endl;
    } else {
    std::cerr<<"Fehler bei dem Erstellen des Streams. Datei vorhanden und lesbar?"<<std::endl;
    return -1;
    }
    return 0;
}
