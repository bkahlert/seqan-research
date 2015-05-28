#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

using namespace std;
using namespace seqan;

int parseLine(char* line){
	int i = strlen(line);
	while (*line < '0' || *line > '9') line++;
		   line[i-3] = '\0';
	i = atoi(line);
	return i;
}


int getValue(){ //Note: this value is in KB!
FILE* file = fopen("/proc/self/status", "r");
int result = -1;
char line[128];


while (fgets(line, 128, file) != NULL){
	if (strncmp(line, "VmSize:", 7) == 0){
		result = parseLine(line);
		break;
	}
}
fclose(file);
return result;
}

int main()
{
	
	
	cout << getValue() << endl;
	String<char> buffer;
	for (int i =0;i<10000000;++i){
		append(buffer,"##########");
	}
	cout << getValue() << endl;
	append(buffer,"###########################");
	cout << getValue() << endl;
	clear(buffer);
	cout << getValue() << endl;
	
	
	
}