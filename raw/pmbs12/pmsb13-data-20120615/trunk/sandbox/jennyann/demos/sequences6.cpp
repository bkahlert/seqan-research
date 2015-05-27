#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

using namespace seqan;
void replaceAs(CharString & myString) {
	typedef Iterator<CharString>::Type it1;
    for (it1 i = begin(myString); i != end(myString); goNext(i))
    {
        if (value(i) == 'a') value(i) = 'X';
        std::cout << value(i);
    }
    std::cout << std::endl;
}

int main()
{
	std::cout << "String: abcdefghijklmnopqrstuvxyz" << std::endl;
	CharString str1 = "abcdefghijklmnopqrstuvxyz";
	replaceAs(str1);
	std::cout << "String: Hello SeqAn!" << std::endl;
	CharString str2 = "Hello SeqAn!";
	replaceAs(str2);
	std::cout << "String: Hello Seqan!" << std::endl;
	CharString str3 = "Hello Seqan!";
	replaceAs(str3);
	return 0;
}