#include <fstream>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/store.h>

using namespace seqan;
int main()
{
	int count = 0; 
	FragmentStore<> store;
	std::ifstream file("assignment_annotations.gtf", std::ios_base::in | std::ios_base::binary);
	read(file, store, Gtf());
	// Create AnnotationTree iterator
	Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
	it = begin(store, AnnotationTree<>());
	// Move iterator one node down
	while (!atEnd(it)&&getType(it) != "exon")goNext(it);
	std::cout << "type: " << getType(it) << endl;
	std::cout << "begin position: " << getAnnotation(it).beginPos << endl;
	std::cout << "end position: " << getAnnotation(it).endPos << endl;
	std::cout << "id: " << value(it) << endl;
	std::cout << "parentid: " << getAnnotation(it).parentId << endl;
	std::cout << "parent name: " << getParentName(it) << endl;
	return 0;
	
}