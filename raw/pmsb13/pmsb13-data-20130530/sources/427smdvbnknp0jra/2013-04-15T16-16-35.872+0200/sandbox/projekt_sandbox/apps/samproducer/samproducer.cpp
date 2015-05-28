#include <fstream>

int main(int argc,char argv)
{	
	std::ofstream myfile;
	myfile.open ("/Informatik/Development/example.sam");
	myfile << "@HD	VN:1.3	SO:coordinate\n@SQ	SN:ref	LN:45\n@SQ	SN:ref2	LN:40\nr001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112\nr002	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*\nr003	0	ref	9	30	5H6M	*	0	0	AGCTAA	*\nr004	0	ref	16	30	6M14N1I5M	*	0	0	ATAGCTCTCAGC	*\nr003	16	ref	29	30	6H5M	*	0	0	TAGGC	*\nr001	83	ref	37	30	9M	=	7	-39	CAGCGCCAT	*";
	myfile.close();

    return 0;
}
