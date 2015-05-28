#include <iostream>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace std;

int main(int argc, char const** args)
{
#pragma omp parallel
	printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());


}
