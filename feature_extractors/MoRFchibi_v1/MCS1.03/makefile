all:
	g++ -O3 source/mcw.c source/*.cpp -o mcw
	g++ -O3 source/mcl.c source/*.cpp -o mcl
	g++ -O3 source/mc.c source/*.cpp -o mc
	g++ source/test.c source/Properties.cpp -o test
	cp inData/Example.fasta input.fasta
clean	:
	rm -f *~ mcw mcl mc test output.txt PSSM/inDir_sp/* PSSM/inDir_ur/* idp/inDir/* MC1/?Results.txt tmp/* source/*~ input.fasta
	rm -rf badSequences
	rm -f clean_input.fasta
	rm -f PSSM/inDir_sp/*
	rm -f PSSM/inDir_ur/*
	rm -f idp/inDir/*
