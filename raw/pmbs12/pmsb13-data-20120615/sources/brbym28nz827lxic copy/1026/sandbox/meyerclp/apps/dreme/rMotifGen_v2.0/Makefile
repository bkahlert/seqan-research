COMPILER = g++
CCFLAGS = -g

rMotifGen: RandomMotifs.cpp SequenceClass.o CommonRoutines.o SubstitutionMatrices.o
	${COMPILER} ${CCFLAGS} -o rMotifGen RandomMotifs.cpp SequenceClass.o CommonRoutines.o SubstitutionMatrices.o

SequenceClass.o: SequenceClass.cpp SequenceClass.H
	${COMPILER} ${CCFLAGS} -c SequenceClass.cpp

SubstitutionMatrices.o: SubstitutionMatrices.cpp SubstitutionMatrices.H
	${COMPILER} ${CCFLAGS} -c SubstitutionMatrices.cpp

CommonRoutines.o: CommonRoutines.H CommonRoutines.cpp SequenceClass.o
	${COMPILER} ${CCFLAGS} -c CommonRoutines.cpp

clean:
	rm *.o a.out rMotifGen