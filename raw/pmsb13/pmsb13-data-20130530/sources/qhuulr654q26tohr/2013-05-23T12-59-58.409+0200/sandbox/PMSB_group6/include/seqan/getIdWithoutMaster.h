//Autor:Hannes

#ifndef GINGER_GETIDWITHOUTMASTER_H_
#define GINGER_GETIDWITHOUTMASTER_H_

#include <seqan/basic.h>

using namespace seqan;
using namespace std;

/**
 * \brief erases the master Id out of a record of a clustered database
 */
CharString getIdWithoutMaster(CharString & id/**<[in]id to be cleaned*/)
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
