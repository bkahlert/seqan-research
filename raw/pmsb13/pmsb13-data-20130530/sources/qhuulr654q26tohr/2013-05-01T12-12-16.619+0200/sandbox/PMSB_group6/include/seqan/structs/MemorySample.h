#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_MEMORYSAMPLE_H
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_MEMORYSAMPLE_H

namespace seqan {
private:
//intern: ??
//prohibited: ??
int parseLine(char* line) {
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}


int getVValue() { //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];


    while (fgets(line, 128, file) != NULL) {
        if (strncmp(line, "VmSize:", 7) == 0) {
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int getPValue() { //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];


    while (fgets(line, 128, file) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}


public:
struct MemorySample {
    CharString name;
    double virtMemory;
    double physMemory;

    MemorySample() {
        computeMemory();
    }
    MemorySample(Charstring n) {
        name=n;
        computeMemory();
    }
    void computeMemory() {
        physMemory=getPValue;
        virtMemory=getVValue;
    }
};
}  // namespace seqan

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_MEMORYSAMPLE_H 
