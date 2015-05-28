
#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/journaled_set.h>

using namespace seqan;
using namespace std;


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
    
    if (length(journalSet) > 0 ) {
        return -1;
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
    
    // [B] Set the reference sequence to the Journal Set
    setGlobalReference(journalSet, tempSeq);
    
    
    // Read remaining sequences.
    while (!atEnd(reader))
    {
        if (readRecord(tempSeqId, tempSeq, reader, Fasta()) != 0)
        {
            std::cerr << "ERROR reading FASTA." << std::endl;
            return 1;
        }

        
        // [C] Append and join the current read sequence.
        appendValue(journalSet, tempSeq);
        join(journalSet, length(journalSet) - 1, joinConfig );
        
    }

    for (int i = 0; i < length(journalSet); ++i) {
        
        cout << value(journalSet, i) << endl;
        

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

}


template <typename TString, typename TPattern>
void findPatternInReference(String<int>& hits, TString const & ref, TPattern const & pattern) {
    
    int last_index = length(ref) - length(pattern) + 1 ;
    for (int i = 0; i < last_index; ++i) {
        bool match = true;
        for (int j = 0; j < length(pattern); ++j) {
            if (ref[i + j] != pattern[j]) {
                match = false;
                break;
            }
        }
        if (match) {
            append(hits, i);       
        }
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
        // Scan the border for a possible match.
        _searchAtBorder(hitTarget, it, journal, pattern);
        ++it;
    }
}     


template <typename TJournalEntriesIterator, typename TPattern>
void _findInOriginalNode(String<int> & hits, TJournalEntriesIterator & it, TPattern const & pattern, String<int> const & refHits) {

    for (int i = 0; i < length(refHits); ++i) {
        int hitPos = refHits[i];

        if (hitPos >= it->physicalPostion && hitPos <= it->physicalPosition + it->length - length(pattern)) {
            int hostPos = hitPos + it->physicalPos - it->virtualPos;
            append(hits, hostPos );            
        } 

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


template <typename TJournalEntriesIterator, typename TPattern>
void _searchAtBorder(String<int> & hits, TJournalEntriesIterator & it, TJournal const & journal,  TPattern const & pattern) {
    
    TJournalIterator patchIter = it->length > lenght(pattern) ? iter(journal, it->virtualPosition + it->length - length(pattern) + 1);
    TJournalIterator pathcEndIter = iter(journal, it->virtualPosition + it->legnth - 1);
    
    while (patchIter != patchEndInter) {
        if (pattern


    }
    

}



int main()
{
    // Definition of the used types.
    typedef String<Dna,Alloc<> > TSequence;
    typedef String<Dna,Journaled<Alloc<>,SortedArray,Alloc<> > > TJournal;
    typedef StringSet< TJournal, Owner<JournaledSet> > TJournaledSet;
    // Open the stream to the file containing the sequences.
    String<char> seqDatabasePath = "sequences.fa";
    std::ifstream databaseFile(toCString(seqDatabasePath), std::ios_base::in);
    if(!databaseFile.good())
    {
        std::cerr << "Cannot open file <" << seqDatabasePath << ">!" << std::endl;
    }
    // Reading each sequence and journal them.
    
    
    // [D] Construct Journaled Set and call loadAndJoin
    TJournaledSet journalSet;
    
    loadAndJoin(journalSet, databaseFile, JoinConfig<GlobalAlign<JournaledCompact> >()  );
    databaseFile.close();
    return 0;
}

