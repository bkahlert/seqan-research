#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/find_motif.h>


using namespace seqan;

int main ()
{
    String<int,int> x;
	addValue(x,((int)1,(int)2));
	std::cout<<"x"<<std::endl;
	return 0;
}