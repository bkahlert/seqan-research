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
    for (unsigned i = 0; i < length(journalSet); ++i)
        findPatternInJournalString(hitSet[i+1], journalSet[i], pattern, hitSet[0]);
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

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPattern>
void findPatternInJournalString(String<int> & hitTarget,
                           String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journal,
                           TPattern const & pattern,
                           String<int> const & refHits)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const TJournal;
    typedef typename JournalType<TJournal>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries>::Type TJournalEntriesIterator;
    if (length(pattern) > length(journal))
        return;
    TJournalEntriesIterator it = begin(journal._journalEntries);
    TJournalEntriesIterator itEnd = findInJournalEntries(journal._journalEntries, length(journal) - length(pattern) + 1) + 1;
    while(it != itEnd)
    {
        if (it->segmentSource == SOURCE_ORIGINAL)
        {   // Find a possible hit in the current source vertex.
            _findInOriginalNode(hitTarget, it, pattern, refHits);
        }
        if (it->segmentSource == SOURCE_PATCH)
        {  // Search for pattern within the patch node.
            _findInPatchNode(hitTarget, it, journal, pattern);
        }
        //Scan the border for a possible match.
        _searchAtBorder(hitTarget, it, journal, pattern);
        ++it;
    }
}

template <typename TJournalEntriesIterator, typename TPattern>
void _findInOriginalNode(String<int> & hitTarget,
                         TJournalEntriesIterator & entriesIt,
                         TPattern const & pattern,
                         String<int> const & refHits)
{
    if (length(refHits)==0)
        return;
    Iterator<String<int> >::Type upper = ::std::upper_bound(begin(refHits), end(refHits), entriesIt->physicalPosition);
    goPrevious(upper);
    if (getValue(upper) < (int)entriesIt->physicalPosition)
        goNext(upper);
    int offset = entriesIt->physicalPosition - entriesIt->virtualPosition;
    while(value(upper) < (int)entriesIt->physicalPosition + (int)entriesIt->length - (int)length(pattern) + 1 && upper != end(refHits))
    {
        appendValue(hitTarget, getValue(upper)+offset);
        goNext(upper);
    }
}

template <typename TJournalEntriesIterator, typename TJournal, typename TPattern>
void _findInPatchNode(String<int> & hitTarget,
                      TJournalEntriesIterator & entriesIt,
                      TJournal const & journal,
                      TPattern const & pattern)
{
    typedef typename Iterator<TJournal const, Standard>::Type TJournalIterator;
// Search for pattern in the insertion node.
    TJournalIterator patchIter = iter(journal, entriesIt->virtualPosition);
    TJournalIterator patchEnd = patchIter + _max(0, (int)entriesIt->length - (int)length(pattern) + 1);
    // Move step by step over search region.
    for (; patchIter != patchEnd; ++patchIter)
    {
        TJournalIterator verifyIter = patchIter;
        bool isHit = true;
        // Search for pattern in the insertion node.
        for (unsigned posPattern = 0; posPattern < length(pattern); ++posPattern, ++verifyIter)
        {
            // Comparing the pattern value with the current value of the iterator.
            if(pattern[posPattern] != getValue(verifyIter))
            {
                isHit = false;
                break;
            }
        }
        if (isHit)
            appendValue(hitTarget, position(patchIter) + entriesIt->virtualPosition);
    }
}

template <typename TJournalEntriesIterator, typename TJournal, typename TPattern>
void _searchAtBorder(String<int> & hitTarget,
                    TJournalEntriesIterator & entriesIt,
                    TJournal const & journal,
                    TPattern const & pattern)
{
    // [A] Determine first position of the at which pattern crosses the border of current node.
    unsigned = entriesIt->physicalPosition + entriesIt->length - length(pattern) + 1;
    // [B] Determine last position before pattern exits the current node or reaches the end of the sequence.
    
    // [C] Move step by step over search region.

    // [D] Scan pattern in current window and report possible hits.
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
    
    String<char> pattern("GTAG");
    
    searchPattern(hitSet, jSet, pattern);
    
    std::cout << globalReference(jSet) << std::endl;
    
    for (unsigned i=0; i<length(hitSet[0]); ++i)
    {
        std::cout << hitSet[0][i] << std::endl;
    }
    return 0;
}