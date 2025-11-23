
//reading from files as cin / cout
const string nmfl="connect";
std::ifstream fileIn (nmfl+".in" );  std::cin.rdbuf(fileIn.rdbuf());
std::ofstream fileOut(nmfl+".out");std::cout.rdbuf(fileOut.rdbuf());

// taking time
auto take_time=[&](){return std::chrono::high_resolution_clock::now();};
auto get_durat=[&](auto start){ return std::chrono::duration_cast<std::chrono::milliseconds>(take_time() - start).count(); };