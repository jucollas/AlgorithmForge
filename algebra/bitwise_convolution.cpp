/*
Autor: Oscar Vargas Pabon

Currently UNTESTED

Note it wont work well for A*A

Taken from https://codeforces.com/blog/entry/119082
*/


template<typename T>
class btw_conv{
public:
	void resize(vector<T> &vec){
		int e=1;while(e<int(vec.size()))e*=2;
		vec.resize(e);
	}

	void trans_subset(vector<T> &vec,int sd){
		// I assume sd\in\{-1,1\}
		for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
			if(i&j)vec[i]+=vec[i^j]*sd;
		}
	}
	void or_conv(vector<T> &A,vector<T> &B){
		// the answer is shown in A; I assume |A|=|B|=2^x for some x
		trans_subset(A,1);trans_subset(B,1);
		rep(i,0,A.size())A[i]*=B[i];
		trans_subset(A,-1);
	}
	
	void trans_superset(vector<T> &vec,int sd){
		// I assume sd\in\{-1,1\}
		for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
			if(i&j)vec[i^j]+=vec[i]*sd;
		}
	}
	
	void and_conv(vector<T> &A,vector<T> &B){
		// the answer is shown in A; I assume |A|=|B|=2^x for some x
		trans_superset(A,1);trans_superset(B,1);
		rep(i,0,A.size())A[i]*=B[i];
		trans_superset(A,-1);
	}
};
