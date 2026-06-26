N ?= 10
opt?=-O2

compile:
	#g++ $(opt) mycmp.cpp -o mfc.out
	g++ $(opt) ../gen.cpp -o gen.out
	g++ $(opt) atc.cpp -o atcoder.out
	g++ $(opt) mont_atc.cpp -o montatc.out
	#g++ $(opt) cp_alg.cpp -o cpalg.out
	#g++ $(opt) recurs.cpp -o recurs.out
	#g++ $(opt) my_iter.cpp -o myiter.out
	g++ $(opt) montg_2rad_clean.cpp -o montgom.out
	#g++ $(opt) ugly_mont.cpp -o uglymont.out
	g++ $(opt) rad4.cpp -o n_rad4.out
	g++ $(opt) weird_recursive.cpp -o weird_rec.out
	g++ $(opt) -DMONTGOMERY weird_recursive.cpp -o weird_rec_montg.out
gen:
	mkdir -p tests
	./gen.out $(N) > tests/test0.in
	./gen.out $(N) > tests/test1.in
	./gen.out $(N) > tests/test2.in
	./gen.out $(N) > tests/test3.in
setup: compile gen
run:
	./run-all.sh
run-all: compile run
do-all: compile gen run
clean:
	rm -r tests
	rm *.out 2> /dev/null #redirects error messages
