echo "--------------------- Start testing flow ------------------"

g++ -std=c++17 -O2 -fsanitize=address,undefined dinitz_test.cpp -o dt.out && ./dt.out && rm dt.out &
g++ -std=c++17 -O2 dinitz_bench.cpp -o db.out && ./db.out --slow && rm db.out &
g++ -std=c++17 -O2 stress_flow_test.cpp -o sf.out && ./sf.out && rm sf.out &
wait
