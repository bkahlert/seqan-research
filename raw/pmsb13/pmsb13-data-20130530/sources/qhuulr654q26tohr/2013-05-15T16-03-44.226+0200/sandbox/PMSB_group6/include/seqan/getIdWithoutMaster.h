#ifndef GINGER_GETIDWITHOUTMASTER_H_
#define GINGER_GETIDWITHOUTMASTER_H_

#include <seqan/basic.h>

using namespace seqan;
using namespace std;

CharString getIdWithoutMaster(CharString & id)
{
	CharString withoutMaster;
	int i=0;
	while (id[i]/='$' && i<length(id)){
		appendValue(withoutMaster, id[i]);
		i++;
	}
	return withoutMaster;
}

#endif  // GINGER_GETIDWITHOUTMASTER_H_
