/*
Autor: Oscar Vargas Pabon

and convolution tested in https://judge.yosupo.jp/problem/bitwise_and_convolution
or convolution currently untested

Note it wont work well for A*A

Taken from https://codeforces.com/blog/entry/119082

Im assuming from my template:
#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
*/


template<typename T>
class btw_conv{
public:
	void resize(vector<T> &vec){
		int e=1;while(e<int(vec.size()))e*=2;
		vec.resize(e);
	}

	void trans_subset(vector<T> &vec,int sd){
		// I assume $sd\in\{-1,1\}$
		for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
			if(i&e)vec[i]+=vec[i^e]*sd;
		}
	}
	void or_conv(vector<T> &A,vector<T> &B){
		// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
		trans_subset(A,1);trans_subset(B,1);
		rep(i,0,A.size())A[i]*=B[i];
		trans_subset(A,-1);
	}
	
	void trans_superset(vector<T> &vec,int sd){
		// I assume $sd\in\{-1,1\}$
		for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
			if(i&e)vec[i^e]+=vec[i]*sd;
		}
	}
	
	void and_conv(vector<T> &A,vector<T> &B){
		// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
		trans_superset(A,1);trans_superset(B,1);
		rep(i,0,A.size())A[i]*=B[i];
		trans_superset(A,-1);
	}
};