#include <iostream>
#include <omp.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace std;

int main(int argc, char const** args)
{
	omp_set_num_threads(4);
#pragma omp parallel
	cout << "Halloah";


}
