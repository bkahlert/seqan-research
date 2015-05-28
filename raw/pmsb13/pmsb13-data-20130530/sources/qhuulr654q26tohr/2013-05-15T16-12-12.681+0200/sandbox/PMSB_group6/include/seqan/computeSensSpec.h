#ifndef GINGER_COMPUTESENSSPEC_H_
#define GINGER_COMPUTESENSSPEC_H_

#include <seqan/basic.h>
#include <seqan/getIdWithoutMaster.h>

using namespace seqan;
using namespace std;

template<typename TScore>
int computeSensSpec(double & sensitivity, double & specitivity, String<GMatch<TScore> > & shouldBe, String<GMatch<TScore> > & actuallyIs, int datalength)
{
	int truePositive=0;
	int trueNegative=0;
	int falsePositive=0;
	int falseNegative=0;
	for (int i=0;i<length(shouldBe);i++){
		for(int j=0;j<length(actuallyIs);j++){
			if(getIdWithoutMaster(shouldBe[i].id)==getIdWithoutMaster(actuallyIs[j].id)&&getIdWithoutMaster(shouldBe[i].queryId)==getIdWithoutMaster(actuallyIs[j].queryId)){
				truePositive++;
			}  
		}
	}
	falseNegative=length(shouldBe)-truePositive;
	falsePositive=length(actuallyIs)-truePositive;
	trueNegative=datalength-falseNegative-falsePositive-truePositive;
	sensitivity=truePositive/(truePositive-falseNegative);
	specitivity=trueNegative/(trueNegative-falsePositive);
	return 0;
}


#endif  // GINGER_COMPUTESENSSPEC_H_
