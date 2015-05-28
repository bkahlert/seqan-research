#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/journaled_set.h>

using namespace seqan;
template <typename TString, typename TStream, typename TSpec>
inline int
loadAndJoin(StringSet<TString, Owner<JournaledSet> > & jSet,
            TStream & stream,
            JoinConfig<TSpec> const & jConfig)
{
    typedef typename Host<TString>::Type THost;
    // Define the RecordReader.
    RecordReader<std::ifstream, SinglePass<> > reader(stream);
    // [A] Ensure the Journal Set is not occupied by other sequences.
    if (length(jSet)>0)
    {
        std::cerr << "Journal set is not empty." << std::endl;
    }
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
    createGlobalReference(jSet, tempSeq);
    // Read remaining sequences.
    while (!atEnd(reader))
    {
        if (readRecord(tempSeqId, tempSeq, reader, Fasta()) != 0)
        {
            std::cerr << "ERROR reading FASTA." << std::endl;
            return 1;
        }
        appendValue(jSet, tempSeq);
    }
    for (unsigned i=0; i<length(jSet); ++i)
    {
        join(jSet, i, jConfig);
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
    String<char> seqDatabasePath = "sequences.fasta";
    std::ifstream databaseFile(toCString(seqDatabasePath), std::ios_base::in);
    if(!databaseFile.good())
    {
        std::cerr << "Cannot open file <" << seqDatabasePath << ">!" << std::endl;
    }
    // Reading each sequence and journal them.
    TJournaledSet jSet;
    
    loadAndJoin(jSet, databaseFile, JoinConfig<GlobalAlign<JournaledCompact> >());

    databaseFile.close();
    return 0;
}