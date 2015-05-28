#include <fstream>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/store.h>

using namespace seqan;
int main()
{ 
	FragmentStore<> store;
	std::ifstream file("assignment_annotations.gtf", std::ios_base::in | std::ios_base::binary);
	read(file, store, Gtf());
	// Create AnnotationTree iterator
	Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
	it = begin(store, AnnotationTree<>());
	// Move iterator one node down
	while (!atEnd(it)&&getType(it) != "exon")
	{
		std::cout << "curr. type: " << getType(it) << std::endl;
		goNext(it);
	}
	std::cout << "type: " << getType(it) <<std::endl ;
	std::cout << "begin position: " << getAnnotation(it).beginPos << std::endl;
	std::cout << "end position: " << getAnnotation(it).endPos << std::endl;
	std::cout << "id: " << value(it) << std::endl;
	std::cout << "parentid: " << getAnnotation(it).parentId <<std:: endl;
	std::cout << "parent name: " << getParentName(it) << std::endl;
	return 0;
	
}