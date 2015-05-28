//Autor:Hannes

#ifndef GINGER_COMPUTESENSSPEC_H_
#define GINGER_COMPUTESENSSPEC_H_

#include <seqan/basic.h>
#include <seqan/getIdWithoutMaster.h>

#include <seqan/GStructs/GMatch.h>
#include <seqan/getIdWithoutMaster.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence>
/**
 * \brief computes sensitivity of the Gsearch algorithm on the clustered database against the nonclustered database.
 * 
 * specifity is can be computed too, therefore the database length is needed, but is commented out of the code since it is not necessary for the GSearch algorithm (Specifity is always 1).
 */
int computeSensSpec(double & sensitivity/**<[out]sensitivity container*/, double & specitivity/**<[out]specifity container*/, String<GMatch<TScore, TSequence> > & shouldBe/**<[in]first result which is to compare*/, String<GMatch<TScore,TSequence> > & actuallyIs/**<[in]second result which is to compare*/, int datalength)
{
	int truePositive=0;
	int trueNegative=0;
	int falsePositive=0;
	int falseNegative=0;
	/*cout << "---------" << endl;
	cout << datalength << endl;
	cout << length(actuallyIs) << endl;
	cout << length(shouldBe) << endl;
	cout << getIdWithoutMaster(shouldBe[0].id) << endl;
	cout << getIdWithoutMaster(shouldBe[1].id) << endl;
	cout << getIdWithoutMaster(shouldBe[2].id) << endl;
	cout << getIdWithoutMaster(actuallyIs[0].id) << endl;
	cout << getIdWithoutMaster(actuallyIs[1].id) << endl;
	cout << getIdWithoutMaster(actuallyIs[2].id) << endl;
	cout << "--" << endl;
	cout << getIdWithoutMaster(shouldBe[0].queryId) << endl;
	cout << getIdWithoutMaster(shouldBe[1].queryId) << endl;
	cout << getIdWithoutMaster(shouldBe[2].queryId) << endl;
	cout << getIdWithoutMaster(actuallyIs[0].queryId) << endl;
	cout << getIdWithoutMaster(actuallyIs[1].queryId) << endl;
	cout << getIdWithoutMaster(actuallyIs[2].queryId) << endl;
	cout << "----" << endl;*/
	for (int i=0;i<length(shouldBe);i++){
		for(int j=0;j<length(actuallyIs);j++){
			//cout << getIdWithoutMaster(shouldBe[i].id) << " == \n" << getIdWithoutMaster(actuallyIs[j].id) << endl;
			//cout << (getIdWithoutMaster(shouldBe[i].id) == getIdWithoutMaster(actuallyIs[j].id)) << endl;
			//cout << length(getIdWithoutMaster(shouldBe[i].id)) << " == \n" << length(getIdWithoutMaster(actuallyIs[j].id)) << endl;
			//cout << getIdWithoutMaster(shouldBe[i].queryId) << " == \n" << getIdWithoutMaster(actuallyIs[j].queryId) << endl;
			//cout << (getIdWithoutMaster(shouldBe[i].queryId) ==  getIdWithoutMaster(actuallyIs[j].queryId)) << endl;
			//cout  << endl;
			if(getIdWithoutMaster(shouldBe[i].id)==getIdWithoutMaster(actuallyIs[j].id)&&getIdWithoutMaster(shouldBe[i].queryId)==getIdWithoutMaster(actuallyIs[j].queryId)){
				truePositive++;
			}  
		}
	}
	falseNegative=length(shouldBe)-truePositive;
	//falsePositive=length(actuallyIs)-truePositive;
	//trueNegative=datalength-falseNegative-falsePositive-truePositive;
	cout << "truePositives:\t" << truePositive << endl;
	cout << "falseNegatives:\t" << falseNegative << endl;
	//cout << "trueNegatives:\t" << trueNegative << endl;
	sensitivity=(double)truePositive/((double)truePositive+falseNegative);
	//specitivity=(double)trueNegative/((double)trueNegative+falsePositive);
	return 0;
}


#endif  // GINGER_COMPUTESENSSPEC_H_
