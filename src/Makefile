../bin/demo_tripod_alt.out: demo_tripod_alt.cc
	mkdir -p ../bin
	g++ -o $@ $< -lgetfem -lboost_system
run:
	make
	../bin/demo_tripod_alt.out
