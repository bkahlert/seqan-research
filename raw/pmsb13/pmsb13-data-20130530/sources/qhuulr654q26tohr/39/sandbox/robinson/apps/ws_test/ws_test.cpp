#include <iostream>
#include <seqan/store.h>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/parallel.h>

using namespace seqan;
using namespace std;


int main(int argc, char const * argv[])
{
	ProfileChar<Dna> pc1;
	ProfileChar<Dna, double> pc2;
	pc1[0]=3;
	pc2[0]=4;
	cout << pc1[0] << endl;
	cout << pc2[0] << endl;
    return 0;
}
