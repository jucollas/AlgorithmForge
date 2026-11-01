 # check the manacher teest
 echo "-----------start running string tests -------------"
 g++ -std=c++17 -O2 -g -D_GLIBCXX_ASSERTIONS -fsanitize=address,undefined manacherTest.cpp -o mt.out && ./mt.out && rm mt.out &
 g++ -std=c++17 -O2 -g -D_GLIBCXX_ASSERTIONS -fsanitize=address,undefined prefFunTest.cpp -o pt.out && ./pt.out && rm pt.out &
 g++ -std=c++17 -O2 -g -D_GLIBCXX_ASSERTIONS -fsanitize=address,undefined zfunTest.cpp -o zt.out && ./zt.out && rm zt.out &
 wait