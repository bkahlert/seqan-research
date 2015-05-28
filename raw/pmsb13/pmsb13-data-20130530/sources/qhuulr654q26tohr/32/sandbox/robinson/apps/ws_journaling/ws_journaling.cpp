#include <iostream>
#include <seqan/file.h>
#include <seqan/journaled_set.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/sequence_journaled.h>

using namespace seqan;
int main() {
	typedef String<Dna, Journaled<Alloc<> , SortedArray, Alloc<> > >
			TJournalString;
	typedef Host<TJournalString>::Type THost;
	typedef StringSet<TJournalString, Owner<JournaledSet> > TJournaledSet;
	TJournaledSet journaledSet;

	std::fstream stream("sequences.fasta", std::ios::binary | std::ios::in);
	if (!stream.good())
		return 1;
	// Read file.
	RecordReader < std::fstream, SinglePass<> > reader(stream);




	CharString id;
	DnaString seq;

	if (!readRecord(id, seq, reader, Fasta()))
		return 1;
	setGlobalReference(journaledSet, seq);

	while (readRecord(id, seq, reader, Fasta())) {
		appendValue(journaledSet, seq);
	}

	for (unsigned i=0; i<length(journaledSet); ++i) {
		join(journaledSet, i, JoinConfig<GlobalAlign<JournaledCompact> > ());
	}
	join(journaledSet, 0, JoinConfig<GlobalAlign<JournaledManhatten> > ()); // Simply inserts the
	join(journaledSet, 1, JoinConfig<GlobalAlign<JournaledCompact> > ()); // Uses default scoring scheme to compute compact journal.

	::std::cout << "Reference: " << globalReference(journaledSet)
			<< ::std::endl;
	for (unsigned i = 0; i < length(journaledSet); ++i)
		::std::cout << "Journaled Sequence " << i << ": " << value(
				journaledSet, i) << ::std::endl;
	return 0;
}
