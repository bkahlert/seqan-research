#include <iostream>
#include <seqan/index.h>
#include <seqan/find.h>

using namespace seqan;


int main()
{
	CharString haystack = "Simon, send more money!";
    StringSet<CharString> needle;
	appendValue(needle, "mo");
	appendValue(needle, "send");
	appendValue(needle, "more");
	
	Finder<CharString> finder(haystack);
    Pattern<StringSet<CharString>, WuManber> pattern(needle);

	Iterator<StringSet<CharString>, Rooted>::Type it = begin(pattern);
	for (; !atEnd(it); goNext(it))
	{
		while (find(finder, value(it)))
        	std::cout << '[' << beginPosition(finder) << ',' << 
			endPosition(finder) << ")\t" << infix(finder) << std::endl;
		//clear()	
	}

	return 0;
}


