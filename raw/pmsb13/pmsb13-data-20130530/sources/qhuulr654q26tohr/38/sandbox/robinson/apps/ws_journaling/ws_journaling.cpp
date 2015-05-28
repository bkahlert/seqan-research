#include <iostream>
#include <seqan/file.h>
#include <seqan/journaled_set.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/sequence_journaled.h>
#include <seqan/find.h>

using namespace std;
using namespace seqan;

typedef String<Dna, Journaled<Alloc<> , SortedArray, Alloc<> > > TJournalString;
typedef Host<TJournalString>::Type THost;
typedef StringSet<TJournalString, Owner<JournaledSet> > TJournaledSet;

template<typename TString, typename TPattern>
void searchPattern(StringSet<String<int> > & hitSet,
		StringSet<TString, Owner<JournaledSet> > const & journalSet,
		TPattern const & pattern) {
	typedef typename Host<TString>::Type THost;
	// Check for valid initial state.
	if (empty(globalReference(journalSet))) {
		::std::cout << "No reference set. Aborted search!" << std::endl;
		return;
	}
	// Reset the hitSet to avoid phantom hits coming from a previous search.
	clear(hitSet);
	resize(hitSet, length(journalSet) + 1);
	// Access the reference sequence.
	THost & globalRef = globalReference(journalSet);
	// Search for pattern in the reference sequence.
	findPatternInSequence(hitSet[0], globalRef, pattern);
	// Search for pattern in the journaled sequences.
	for (unsigned i = 0; i < length(journalSet); ++i)
		findPatternInSequence(hitSet[i + 1], journalSet[i], pattern,
				hitSet[0]);
}

template<typename TString, typename TPattern>
void findPatternInSequence(String<int> & hits, TString const & referenceSeq,
		TPattern const & patternSeq) {
	Finder<TString> finder(referenceSeq);
	Pattern<TPattern, Horspool> pattern(patternSeq);
	while (find(finder, pattern))
		++hits[beginPosition(finder)];
}

int main() {

	TJournaledSet journaledSet;

	std::fstream stream("sequences.fasta", std::ios::binary | std::ios::in);
	if (!stream.good()) {
		cerr << "could not open file" << endl;
		return 1;
	}
	// Read file.
	RecordReader < std::fstream, SinglePass<> > reader(stream);

	CharString id;
	DnaString seq;

	if (readRecord(id, seq, reader, Fasta())) {
		cerr << "could not read reference" << endl;
		return 1;
	}
	setGlobalReference(journaledSet, seq);

	while (!readRecord(id, seq, reader, Fasta())) {
		appendValue(journaledSet, seq);
	}

	for (unsigned i = 0; i < length(journaledSet); ++i) {
		join(journaledSet, i, JoinConfig<GlobalAlign<JournaledCompact> > ());
	}

	::std::cout << "Reference: " << globalReference(journaledSet)
			<< ::std::endl;
	for (unsigned i = 0; i < length(journaledSet); ++i)
		::std::cout << "Journaled Sequence " << i << ": " << value(
				journaledSet, i) << ::std::endl;
	return 0;
}
