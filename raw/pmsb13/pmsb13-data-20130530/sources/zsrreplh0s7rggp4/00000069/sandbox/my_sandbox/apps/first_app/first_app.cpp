#include <seqan/sequence.h>
#include <seqan/index.h>
#include <iostream>
#include <seqan/sequence_journaled.h>
#include <seqan/journaled_set.h>

using namespace std;
using namespace seqan;


template <typename TString, typename TStream, typename TSpec>
inline int
loadAndJoin(StringSet<TString, Owner<JournaledSet> > & journalSet,
            TStream & stream,
            JoinConfig<TSpec> const & joinConfig)
{
    typedef typename Host<TString>::Type THost;
    // Define the RecordReader.
    RecordReader<std::ifstream, SinglePass<> > reader(stream);
    // [A] Ensure the Journal Set is not occupied by other sequences.
    clear(journalSet);
    // Construct the temporary buffers for the read id and sequence.
    String<char> tempSeqId;
    THost tempSeq;
    // No sequences in the fasta file!
    if (atEnd(reader))
    {
        std::cerr << "Empty FASTA file." << std::endl;
        return -1;
    }
    // First read sequence for reference sequence.
    if (readRecord(tempSeqId, tempSeq, reader, Fasta()) != 0)
    {
        std::cerr << "ERROR reading FASTA." << std::endl;
        return 1;
    }
    // [B] Set the reference sequence to the Journal Set
    createGlobalReference(journalSet, tempSeq);
    // Read remaining sequences.
    while (!atEnd(reader))
    {
        if (readRecord(tempSeqId, tempSeq, reader, Fasta()) != 0)
        {
            std::cerr << "ERROR reading FASTA." << std::endl;
            return 1;
        }
        // [C] Append and join the current read sequence.
    appendValue(journalSet,tempSeq);
    join(journalSet,length(journalSet)-1,joinConfig);  
      
    }
    return 0;
}
int main()
{
    // Definition of the used types.
    typedef String<Dna,Alloc<> > TSequence;
    typedef String<Dna,Journaled<Alloc<>,SortedArray,Alloc<> > > TJournal;
    typedef StringSet< TJournal, Owner<JournaledSet> > TJournaledSet;
    // Open the stream to the file containing the sequences.
    String<char> seqDatabasePath = "/home/gabriel/Downloads/sequences.fasta";
    std::ifstream databaseFile(toCString(seqDatabasePath), std::ios_base::in);
    if(!databaseFile.good())
    {
        std::cerr << "Cannot open file <" << seqDatabasePath << ">!" << std::endl;
    }
    // Reading each sequence and journal them.
    // [D] Construct Journaled Set and call loadAndJoin
    TJournaledSet journalSet;
    JoinConfig<GlobalAlign<JournaledCompact> > joinConfig;
    loadAndJoin(journalSet, databaseFile, joinConfig);
    databaseFile.close();
    return 0;
}