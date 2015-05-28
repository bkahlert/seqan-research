//Autor:Jakob

#ifndef GINGER_COMPAREFILES_H_
#define GINGER_COMPAREFILES_H_

#include <seqan/basic.h>

using namespace seqan;
using namespace std;

/** \brief returns 0 if two files are exactly the same*/
int compareFiles(FILE* f1, FILE* f2) {
	int N = 10000;
	char buf1[N];
	char buf2[N];
	
	do {
		size_t r1 = fread(buf1, 1, N, f1);
		size_t r2 = fread(buf2, 1, N, f2);
		
		if (r1 != r2 ||
			memcmp(buf1, buf2, r1)) {
			SEQAN_FAIL("Content is not the same."); // Files are not equal
			return 1;
			}
	} while (!feof(f1) && !feof(f2));
	if (!(feof(f1) && feof(f2))){
		SEQAN_FAIL("File length is not the same.");
		return 1;
	}
	return 0;
}

#endif  // GINGER_COMPAREFILES_H_
