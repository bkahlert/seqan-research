cmake_minimum_required (VERSION 2.6)
project (core_apps_rabema)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(rabema rabema.cpp
			    rabema.h
			    build_gold_standard.h
			    build_gold_standard_options.h
			    curve_smoothing.h
			    evaluation.h
			    evaluation_options.h
			    find_approx_dp_quality.h
			    find_hamming_simple_ext.h
			    find_hamming_simple_quality.h
			    find_myers_ukkonen_ext.h
			    find_myers_ukkonen_reads.h
			    intervals.h
			    verification.h
			    wit_store.h
			    witio.h)

seqan_add_executable(do_search do_search.cpp witio.h find_myers_ukkonen_ext.h)

