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
    clear(jSet);
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

template <typename TString, typename TPattern>
void searchPattern(StringSet<String<int> > & hitSet,
                   StringSet<TString, Owner<JournaledSet> > const & journalSet,
                   TPattern const & pattern)
{
    typedef typename Host<TString>::Type THost;
    // Check for valid initial state.
    if (empty(globalReference(journalSet)))
    {
        ::std::cout << "No reference set. Aborted search!" << std::endl;
        return;
    }
    // Reset the hitSet to avoid phantom hits coming from a previous search.
    clear(hitSet);
    resize(hitSet, length(journalSet) + 1);
    // Access the reference sequence.
    THost & globalRef = globalReference(journalSet);
    // Search for pattern in the reference sequence.
    findPatternInReference(hitSet[0], globalRef, pattern);
    // Search for pattern in the journaled sequences.
//    for (unsigned i = 0; i < length(journalSet); ++i)
//        findPatternInJournalString(hitSet[i+1], journalSet[i], pattern, hitSet[0]);
}

template <typename TString, typename TPattern>
void findPatternInReference(String<int> & hits,
                            TString const & reference,
                            TPattern const & pattern)
{
    if (length(pattern) > length(reference))
    {
        return;
    }
    for (unsigned i=0; i<length(reference)-length(pattern); ++i)
    {
        unsigned j=0; 
        while (j<length(pattern) && reference[i+j] == pattern[j])
            ++j; 
        if (j==length(pattern))
            appendValue(hits, i);
    }

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
    
    StringSet<String<int> > hitSet;
    
    String<char> pattern("AAAAAA");
    
    searchPattern(hitSet, jSet, pattern);
    
    return 0;
}