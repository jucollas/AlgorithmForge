@echo off

g++ -O2 mycmp.cpp -o mfc
g++ -O2 ..\gen.cpp -o gen

g++ -O2 atc.cpp -o atcoder
g++ -O2 mont_atc.cpp -o montatc
::g++ -O2 cp_alg.cpp -o cpalg
::g++ -O2 recurs.cpp -o recurs
::g++ -O2 my_iter.cpp -o myiter
g++ -O2 montg_2rad_clean.cpp -o montgom
::g++ -O2 ugly_mont.cpp -o uglymont

g++ -O2 rad4.cpp -o n_rad4

echo on