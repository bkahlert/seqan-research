//Autor:Hannes

#ifndef GINGER_GETIDWITHOUTMASTER_H_
#define GINGER_GETIDWITHOUTMASTER_H_

#include <seqan/basic.h>

using namespace seqan;
using namespace std;

CharString getIdWithoutMaster(CharString & id)
{
	CharString withoutMaster;
	int i=0;
	while (i<length(id) && id[i]!='$'){
		appendValue(withoutMaster, id[i]);
		i++;
	}
	return withoutMaster;
}

#endif  // GINGER_GETIDWITHOUTMASTER_H_
