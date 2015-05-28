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
	};
}  // namespace seqan

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_PERFORMANCESAMPLE_H 
