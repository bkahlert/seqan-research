#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>


using namespace seqan;

int main ()
{
    String<int,int> x;
	addeValue(x,<1,2>);
	std::cout<<"x"<<endl;
	return 0;
}