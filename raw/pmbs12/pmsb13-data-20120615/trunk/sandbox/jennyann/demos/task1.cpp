#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>

using namespace seqan;

template <typename TAlphabet>
void countchar(TAlphabet const &)
{
	String<char> my_AminoAcids = "MQDRVKRPMNAFIVWSRDQRRKMALEN";

	typedef Iterator<String<char>>::Type it1;
    for (it1 i = begin(my_AminoAcids); i != end(my_AminoAcids); goNext(i))
    {
        if (value(i) == 'R') value(i) = 'A';
        std::cout << value(i) << ',';
    }
    std::cout << std::endl;

	typedef typename Size<String<char>>::Type TStringSize;
    TStringSize strLen = length(my_AminoAcids);
    typedef typename Size<TAlphabet>::Type TSize;
    TSize alphSize = ValueSize<TAlphabet>::VALUE;

	int counter = 0;

    for (TSize i = 0; i < alphSize; ++i)
	{
		for (TStringSize j = 0; j < strLen; ++j)
		{
			if (value(my_AminoAcids, j) == TAlphabet(i))
			{
				counter += 1;
			}
		}
		std::cout << TAlphabet(i) << ',' << counter << '\n';
		counter = 0;
	}
	
    std::cout << std::endl;
}

int main()
{
	countchar(AminoAcid());
	return 0;
}