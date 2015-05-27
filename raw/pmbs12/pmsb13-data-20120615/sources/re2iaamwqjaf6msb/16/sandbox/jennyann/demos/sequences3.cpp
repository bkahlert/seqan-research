#include <iostream>

#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

void replaceChar(String<char> str, int pos)
{
	for (char c = 'a'; c <= 'z'; ++c) {
		str[pos] = c;
		std::cout << str << std::endl;
	}
}

void printPermutations(int len)
{
	String<char> str;
  for (int i = 0; i < len; ++i)
  {
    str += "a";
  }
  std::cout << str << std::endl;

  for (int i = 0; i < len; ++i)
  {
	  replaceChar(str, i);
  }

}

int main() {
  printPermutations(3);
  return 0;
}