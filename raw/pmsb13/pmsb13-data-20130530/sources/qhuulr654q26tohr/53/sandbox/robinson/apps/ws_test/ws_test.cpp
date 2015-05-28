#include <iostream>
#include <fstream>
#include <seqan/store.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/parallel.h>
#include </home/Development/seqan-trunk/sandbox/robinson/include/seqan/readJasper.h>


using namespace seqan;
using namespace std;


int main(int argc, char const ** argv) {
	// Handle command line arguments, open files.
	if (argc != 2)
		return 1;
	std::fstream stream(argv[1], std::ios::binary | std::ios::in);
	if (!stream.good())
		return 1;

	// Read file.
	RecordReader < fstream, SinglePass<> > reader(stream);
	CharString name;
	CharString id;
	String < ProfileChar<Dna> > matrix;

	int res = readRecord(id, name, matrix, reader, Jaspar());
	if (res != 0)
		return res;

	// Write out some of the data to stdout.
	cout << id << "\t" << name << endl;
	cout << matrix[0].count[0] << endl;
	cout << matrix[0].count[1] << endl;
	cout << matrix[0].count[2] << endl;
	cout << matrix[0].count[3] << endl;
	return 0;
}
