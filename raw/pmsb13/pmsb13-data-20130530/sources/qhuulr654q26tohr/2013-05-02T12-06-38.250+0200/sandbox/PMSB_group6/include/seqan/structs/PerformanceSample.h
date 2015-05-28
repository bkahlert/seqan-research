#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_PERFORMANCESAMPLE_H
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_PERFORMANCESAMPLE_H

namespace seqan {
	struct PerformanceSample {
		CharString name;
		double startTime;
		double endTime;
		
		PerformanceSample() {
			start();
		}
		PerformanceSample(CharString n) :
		name(n) {
			start();
		}
		PerformanceSample(CharString n, double start) :
		name(n), startTime(start) {
		}
		void start(){
			startTime=sysTime();
		}
		void end(){
			endTime=sysTime();
		}
		double getTime(){
			return (endTime-startTime);
		}
		void printToStdout() {
			std::cout << "Runtime - " << name << std::endl;
			std::cout << getTime() << " s" <<  std::endl << std::endl;
		}
		void takeAverageValues(String<PerformanceSample> list){
			startTime=0.0;
			endTime=0.0;
			for (int i=0;i<length(list);++i){
				startTime+=list[i].startTime;
				endTime+=list[i].endTime;
			}
			endTime=endTime-startTime;
			startTime=0;
		}
	};
}  // namespace seqan

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_PERFORMANCESAMPLE_H 
