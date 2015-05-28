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

template <typename TString, typename TPattern>
void findPatternInReference(String<int> & hits,
                            TString const & reference,
                            TPattern const & pattern)
{
    // [A] Check whether pattern fits into the sequence.
    if (length(pattern) > length(reference))
        return;
    // [B] Iterate over all positions at which the pattern might occur.
    for (unsigned pos = 0; pos < length(reference) - length(pattern) + 1; ++pos)
    {
        bool isHit = true;
        // [C] Evaluate all positions of the pattern until you find a mismatch or you have found a hit.
        for (unsigned posPattern = 0; posPattern < length(pattern); ++posPattern)
        {
            if(pattern[posPattern] != reference[posPattern + pos])
            {
                isHit = false;
                break;
            }
        }
        // [D] Report begin position at which pattern matches the sequence.
        if(isHit)
            appendValue(hits, pos);
    }
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