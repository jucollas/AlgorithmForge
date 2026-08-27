/* Author: Oscar Vargas Pabon
based loosely on https://atcoder.jp/contests/arc222/submissions/76586635

Tested locally. Remember to do the /n in the end for XOR.
For 'and','or', sd=1 means forward transformation and sd=0 means inverse tranformation
*/ template <typename T, int n> void fwt_rc(T*A){
	if constexpr(n==1) return;
	constexpr int m=n/2;
	fwt_rc<T,m>(A);fwt_rc<T,m>(A+m);
	for(int i=0;i<m;++i){
		const T x=A[i],y=A[i+m];
		A[i]=x+y;A[i+m]=x-y;
		// A[i] = (sd)?x+y:x-y;    // and
		// A[i+m] = (sd)?x+y:y-x; // or
	}
}template<typename T,int n=1<<25> void fwt(vector<T>&A){
	if constexpr(n==1)return;//.data() seems powerfull
	if(int(A.size())==n) fwt_rc<T,n>(A.data());
	else fwt<T,n/2>(A);
}template<typename T> void convolution(vector<T> &A, vector<T> B) {
	const int n=A.size(); assert(A.size()==B.size());
	fwt(A); fwt(B);
	for(int i=0;i<n;++i)A[i]*=B[i];
	fwt(A); const T tn=T(n).inv();
	for(int i=0;i<n;++i)A[i]*=tn; }