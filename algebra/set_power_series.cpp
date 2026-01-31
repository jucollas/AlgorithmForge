/*
Author: Oscar Vargas Pabon

and convolution tested in https://judge.yosupo.jp/problem/bitwise_and_convolution
or convolution currently untested

exp and log tested in https://judge.yosupo.jp/submission/349366

All impl are also tested in my personal libs

Taken from https://codeforces.com/blog/entry/119082
			https://codeforces.com/blog/entry/92128
			https://gist.github.com/dario2994/e3257326ee80c054d3b48766b600991a

Im assuming from my template:
#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
*/

template<typename T>
void trans_subset(vector<T> &vec,int sd){
	// I assume $sd\in\{-1,1\}$
	// When sd=+1 $vec[k]'=\sum_{i,j\subseteq k}vec[i] *vec[j] $
	// When sd=-1 $vec[k] =\sum_{i,j\subseteq k}vec[i]'*vec[j]'$
	for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
		if(i&e)vec[i]+=vec[i^e]*T(sd);
	}
}
template<typename T>
vector<T> or_conv(vector<T> A,vector<T> B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes $A'[k]=\sum_{(i|j)==k}A[i]*B[j]$
	trans_subset(A,1);trans_subset(B,1);
	rep(i,0,A.size())A[i]*=B[i];
	trans_subset(A,-1);
	return A;
}
template<typename T>	
void trans_superset(vector<T> &vec,int sd){
	// I assume $sd\in\{-1,1\}$
	// When sd=+1 $vec[k]'=\sum_{i,j\superseteq k}vec[i] *vec[j] $
	// When sd=-1 $vec[k] =\sum_{i,j\superseteq k}vec[i]'*vec[j]'$
	for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
		if(i&e)vec[i^e]+=vec[i]*T(sd);
	}
}
template<typename T>	
vector<T> and_conv(vector<T> &A,vector<T> &B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes A'[k]=\sum_{(i&j)==k}A[i]*B[j]$
	trans_superset(A,1);trans_superset(B,1);
	rep(i,0,A.size())A[i]*=B[i];
	trans_superset(A,-1);
	return A;
}
template<typename T>
vector<T> subset_conv(const vector<T> &A,const vector<T>&B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes $C[k]=\sum_{s\subseteq k}A[s]*B[k\setminus s]$
	const int n=ilog2(A.size()); assert(A.size()==B.size()&&(1<<n)==int(A.size()));
	vector<vector<T>> A_hat(n+1,vector<T>(1<<n,0)),B_hat=A_hat,C_hat=B_hat;
	rep(i,0,1<<n){ int pcnt=__builtin_popcount(i);
		A_hat[pcnt][i]=A[i]; B_hat[pcnt][i]=B[i];
	} rep(i,0,n+1)trans_subset(A_hat[i],1),trans_subset(B_hat[i],1);
	rep(k,0,n+1)rep(i,0,k+1)rep(j,0,1<<n){
		C_hat[k][j]+=A_hat[i][j]*B_hat[k-i][j];
	}rep(i,0,n+1)trans_subset(C_hat[i],-1);
	vector<T> C(1<<n);rep(i,0,1<<n){
		C[i]=C_hat[__builtin_popcount(i)][i];
	}return C;
}
template<typename T> vector<T> subset_iconv(const vector<T>&A,const vector<T>&C){
	// Computes B such that $C[k]=\sum_{i\subseteq k}A[i]*B[k\setminus i]$
	const T A0_inv=A[0].inv(); assert(!(A[0]==0));// Im assuming some .inv()
	const int n=ilog2(A.size()); assert(A.size()==C.size()&&(1<<n)==int(A.size()));
	vector<vector<T>> A_hat(n+1,vector<T>(1<<n,0)),B_hat=A_hat;
	rep(i,0,1<<n)A_hat[__builtin_popcount(i)][i]=A[i];
	rep(i,0,n+1)trans_subset(A_hat[i],1);
	vector<T>B(1<<n); rep(k,0,n+1){
		rep(i,0,k)rep(j,0,1<<n)B_hat[k][j]+=B_hat[i][j]*A_hat[k-i][j];
		trans_subset(B_hat[k],-1);
		rep(j,0,1<<n){
			// C_hat[k][j] = B_hat[k] =0 since ppcnt(j)!=k
			if(__builtin_popcount(j)!=k)B_hat[k][j]=0;
			//C_hat[k][j]=C[j]=B[j]*A[0]+B_hat[k][j]
			else B[j]=B_hat[k][j]=(C[j]-B_hat[k][j])*A0_inv;
		} trans_subset(B_hat[k],1);
	} return B;
}
template<typename T> vector<T> set_exp(const vector<T>&A){
	// Computes $exp(A)=\sum_{i\geq 0}A^i/i!$
	// exp(A)[k]=\sum_{p\in P[k]}\prod_{i\in p}A_{p_i}$
	// Where P is the set of all partitions of k into nonempty sets
	assert(A[0]==0);// required for A^i to be nilpotent
	vector<T> B={1};B.reserve(A.size());
	for(int e=1;e<int(A.size());e<<=1){
		for(const T &ac:subset_conv({A.begin()+e,A.begin()+e*2},B)){
			B.emplace_back(ac);
		}
	}return B;//B=exp(A)
}
template<typename T> vector<T> set_log(const vector<T>&A){
	// Computes given A=exp(B) ; computes B;
	assert(A[0]==1);
	vector<T>B={0};B.reserve(A.size());
	for(int e=1;e<int(A.size());e<<=1){
		for(const T &ac:subset_iconv<T>( {A.begin(),A.begin()+e},{A.begin()+e,A.begin()+e*2} )){
			B.emplace_back(ac);
		}
	}return B;
}