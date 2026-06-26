/*
Shamelessly taken from https://github.com/mochow13/competitive-programming-library/blob/master/Math/Fast%20Walsh-Hadamard%20Transform.cpp
and adapted by yours truly, osvarp

Tested in https://codeforces.com/contest/1411/problem/G
*/

// untested
// based loosely on https://atcoder.jp/contests/arc222/submissions/76586635
template <typename T, int n>
void fwt_rc(T*A){
	if constexpr(n==1) return;
	constexpr int m=n/2;
	fwt_rc<T,m>(A);fwt_rc<T,m>(A+m);
	for(int i=0;i<n;++i){
		const T&x=A[i],&y=A[i+m];
		A[i]=x+y;A[i+m]=x-y;
		// A[i] = (sd)?x+y:x-y;    // and
		// A[i+m] = (sd)?x+y:y-x; // or
	}
}template<typename T,int n=1<<20>
void fwt(vector<T>&A){
	if constexpr(n==1)return;//.data() seems powerfull
	if(int(A.size())==n)return fwt_rc<T,n>(A.data());
	return fwt<T,n/2>(A);
} void convolution(vector<T> &A, vector<T> B) {
	const int n=A.size(); assert(A.size()==B.size());
	fwt(A); fwt(B);
	for(int i=0;i<n;++i)A[i]*=B[i];
	fwt(A);
	const T tn=n;
	for(int i=0;i<n;++i)A[i]/=tn;
}

template <typename T>
struct FWT {
	void fwt( vector<T> & io, bool sd ) {
		const int n=io.size();
		for (int d = 1; d < n; d <<= 1) {
			for (int i = 0, m = d<<1; i < n; i += m) {
				for (int j = 0; j < d; j++) { /// Don't forget modulo if required
					T x = io[i+j], y = io[i+j+d];
					
					io[i+j] = (x+y), io[i+j+d] = (x-y);	// xor
					if(!sd)io[i+j]/=2,io[i+j+d]/=2;/// Modular inverse if required here
					
					// io[i+j] = (sd)?x+y:x-y;    // and
					// io[i+j+d] = (sd)?x+y:y-x; // or
				}
			}
		}
	}
	// a, b are two polynomials and n is size which is power of two
	void convolution(vector<T> &A, vector<T> &B) {
		const int n=A.size();
		fwt(A,1); fwt(B,1);
		for (int i = 0; i < n; i++)
			A[i] = A[i]*B[i];/// Don't forget modulo if required
		fwt(A,0);
	}
	// for a*a	
	void self_convolution(vector<T> &A) {
		const int n=A.size();
		fwt(A,1);
		for (int i = 0; i < n; i++)
			A[i] = A[i]*A[i];/// Don't forget modulo if required
		fwt(A,0);
	}
};FWT<mint> fwt;