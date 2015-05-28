#include <iostream>
#include <seqan/file.h>
#include <seqan/journaled_set.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/sequence_journaled.h>


using namespace std;
using namespace seqan;

int main() {
	typedef String<Dna, Journaled<Alloc<> , SortedArray, Alloc<> > >
			TJournalString;
	typedef Host<TJournalString>::Type THost;
	typedef StringSet<TJournalString, Owner<JournaledSet> > TJournaledSet;
	TJournaledSet journaledSet;

	std::fstream stream("sequences.fasta", std::ios::binary | std::ios::in);
	if (!stream.good()){
		cerr << "could not open file" << endl;
		return 1;}
	// Read file.
	RecordReader < std::fstream, SinglePass<> > reader(stream);




	CharString id;
	DnaString seq;

	if (readRecord(id, seq, reader, Fasta())){
		cerr << "could not read reference" << endl;
		return 1;}
	setGlobalReference(journaledSet, seq);

	while (readRecord(id, seq, reader, Fasta())) {
		appendValue(journaledSet, seq);
	}

	for (unsigned i=0; i<length(journaledSet); ++i) {
		join(journaledSet, i, JoinConfig<GlobalAlign<JournaledCompact> > ());
	}

	::std::cout << "Reference: " << globalReference(journaledSet)
			<< ::std::endl;
	for (unsigned i = 0; i < length(journaledSet); ++i)
		::std::cout << "Journaled Sequence " << i << ": " << value(
				journaledSet, i) << ::std::endl;
	return 0;
}
