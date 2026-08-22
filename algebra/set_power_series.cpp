/* Author: Oscar Vargas Pabon
and convolution tested in https://judge.yosupo.jp/problem/bitwise_and_convolution
or convolution currently untested
exp and log tested in https://judge.yosupo.jp/submission/396113
All impl are also tested in my personal libs
Taken from https://codeforces.com/blog/entry/119082
			https://codeforces.com/blog/entry/92128
			https://gist.github.com/dario2994/e3257326ee80c054d3b48766b600991a
Note the not-subnormal way to do this is (the first is subset, second superset)
for(int e=1;e<n;e*=2)for(int i=0;i<n;++i)if(i&e)vec[i  ]+=tp?vec[i^e]:-vec[i^e];
for(int e=1;e<n;e*=2)for(int i=0;i<n;++i)if(i&e)vec[i^e]+=tp?vec[i  ]:-vec[i  ];
REMEMBER CHECKING max_exp !!!
*/ const int max_exp=22;// (1<<max_exp)> maximum array size possible
template<typename tint>void raw_subset(tint*vec,int n,bool tp){
	if(tp)
	for(int e=1;e<n;e*=2)for(int r=0;r<n;r+=2*e)for(int l=0;l<e;++l)
		vec[l|r|e]+=vec[l|r];
	else 
	for(int e=1;e<n;e*=2)for(int r=0;r<n;r+=2*e)for(int l=0;l<e;++l)
		vec[l|r|e]-=vec[l|r];
}template<typename tint>inline void subset(vector<tint>&vec,bool tp){
	// When tp=1 $vec[k]'=\sum_{i,j\subseteq k}vec[i] *vec[j] $
	// When tp=0 $vec[k] =\sum_{i,j\subseteq k}vec[i]'*vec[j]'$
	raw_subset<tint>(vec.data(),vec.size(),tp);
}template<typename tint>void raw_superset(tint*vec,int n,bool tp){
	if(tp)
	for(int e=1;e<n;e*=2)for(int r=0;r<n;r+=2*e)for(int l=0;l<e;++l)
		vec[l|r]+=vec[l|r|e];
	else
	for(int e=1;e<n;e*=2)for(int r=0;r<n;r+=2*e)for(int l=0;l<e;++l)
		vec[l|r]-=vec[l|r|e];
}template<typename tint>inline void superset(vector<tint>&vec,bool tp){
	// When tp=1 $vec[k]'=\sum_{i,j\superseteq k}vec[i] *vec[j] $
	// When tp=0 $vec[k] =\sum_{i,j\superseteq k}vec[i]'*vec[j]'$
	raw_superset<tint>(vec.data(),vec.size(),tp);
} template<typename tint>vector<tint>or_conv(vector<tint>A,vector<tint>B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes $A'[k]=\sum_{(i|j)==k}A[i]*B[j]$
	subset<tint>(A,1);subset<tint>(B,1);
	for(int i=0;i<int(A.size());++i)A[i]*=B[i];
	subset<tint>(A,0); return A;
} template<typename tint>vector<tint>and_conv(vector<tint>A,vector<tint>B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes A'[k]=\sum_{(i&j)==k}A[i]*B[j]$
	superset<tint>(A,1);superset<tint>(B,1);
	for(int i=0;i<int(A.size());++i)A[i]*=B[i];
	superset<tint>(A,0); return A;
} template<typename tint>void raw_subset_conv(const tint A[],const tint B[],tint*C,int n){
	// I assume $2^n=|A|=|B|=|C|$
	static tint A_hat[(max_exp+1)<<max_exp],B_hat[(max_exp+1)<<max_exp],C_hat[(max_exp+1)<<max_exp];
	for(int i=0;i<((n+1)<<n);++i)A_hat[i]=B_hat[i]=C_hat[i]=0;
	for(int i=0;i<(1<<n);++i){int pcnt=__builtin_popcount(i);
		A_hat[pcnt<<n|i]=A[i]; B_hat[pcnt<<n|i]=B[i];
	}for(int i=0;i<=n;++i)raw_subset<tint>(A_hat+(i<<n),1<<n,1),
	                      raw_subset<tint>(B_hat+(i<<n),1<<n,1);
	for(int k=0;k<=n;++k)for(int i=0;i<=k;++i)for(int j=0;j<(1<<n);++j)
		C_hat[k<<n|j]+=A_hat[i<<n|j]*B_hat[(k-i)<<n|j];
	for(int i=0;i<=n;++i)raw_subset<tint>(C_hat+(i<<n),1<<n,0);
	for(int i=0;i<(1<<n);++i)C[i]=C_hat[__builtin_popcount(i)<<n|i];
}template<typename tint>inline vector<tint>subset_conv(const vector<tint>&A,const vector<tint>&B){
	// Computes $C[k]=\sum_{s\subseteq k}A[s]*B[k\setminus s]$
	const int n=ilog2(A.size()); assert(A.size()==B.size()&&(1<<n)==int(A.size()));
	vector<tint> C(1<<n); raw_subset_conv<tint>(A.data(),B.data(),C.data(),n);
	return C;
} template<typename tint>void raw_subset_iconv(const tint A[],const tint C[],tint B[],int n){
	// I assume $2^n=|A|=|B|=|C|$
	const tint A0_inv=A[0].inv(); assert(!(A[0]==0));// Im assuming some .inv()
	static tint A_hat[(max_exp+1)<<max_exp],B_hat[(max_exp+1)<<max_exp];
	for(int i=0;i<((n+1)<<n);++i)A_hat[i]=B_hat[i]=0;
	for(int i=0;i<(1<<n);++i)A_hat[__builtin_popcount(i)<<n|i]=A[i];
	for(int i=0;i<=n;++i)raw_subset<tint>(A_hat+(i<<n),1<<n,1);
	for(int k=0;k<=n;++k){
		for(int i=0;i<k;++i)for(int j=0;j<(1<<n);++j)
			B_hat[k<<n|j]+=B_hat[i<<n|j]*A_hat[(k-i)<<n|j];
		raw_subset<tint>(B_hat+(k<<n),1<<n,0);
		for(int j=0;j<(1<<n);++j){
			if(__builtin_popcount(j)!=k)B_hat[k<<n|j]=0;
			else B[j]=B_hat[k<<n|j]=(C[j]-B_hat[k<<n|j])*A0_inv;
		} raw_subset<tint>(B_hat+(k<<n),1<<n,1); }
}template<typename tint>inline vector<tint>subset_iconv(const vector<tint>&A,const vector<tint>&C){
	// Computes B such that $C[k]=\sum_{i\subseteq k}A[i]*B[k\setminus i]$
	const int n=ilog2(A.size()); assert(A.size()==C.size()&&(1<<n)==int(A.size()));
	vector<tint>B(1<<n);raw_subset_iconv<tint>(A.data(),C.data(),B.data(),n);return B;
} template<typename tint>vector<tint>set_exp(const vector<tint>&A){
	// Computes $exp(A)=\sum_{i\geq 0}A^i/i!$
	// exp(A)[k]=\sum_{p\in P[k]}\prod_{i\in p}A_{p_i}$
	// Where P is the set of all partitions of k into nonempty sets
	assert(A[0]==0);// required for A^i to be nilpotent
	vector<tint> B(A.size());B[0]=1; for(int e=0;(1<<e)<int(A.size());++e)
		raw_subset_conv<tint>(A.data()+(1<<e),B.data(),B.data()+(1<<e),e);
	return B;//B=exp(A)
}template<typename tint>vector<tint>set_log(const vector<tint>&A){
	assert(A[0]==1); // Computes given A=exp(B) ; computes B;
	vector<tint>B(A.size(),0); for(int e=0;(1<<e)<int(A.size());++e)
		raw_subset_iconv<tint>(A.data(),A.data()+(1<<e),B.data()+(1<<e),e);
	return B; }