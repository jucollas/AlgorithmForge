/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe en testing\rangeQtest.cpp
*/
template<typename tint> struct fenw{
	vector<tint>bit; int n;
	fenw(int nn=0){ n=nn;bit.resize(n+1,0);}
	fenw(const vector<tint>&arr) { assert(arr[0]==0);//note arr[0] is undefined
		n=arr.size();bit.resize(n+1,0);fill(bit.begin(),bit.end(),tint(0));
		int nxt; for(int i=0;i<=n;++i) {
			bit[i]+=arr[i];nxt=i+(i&(-i));
			if(nxt<=n)bit[nxt]+=bit[i]; }
	} void update( int x, tint val ) {
		for(;x<=n;x+=x&(-x))bit[x]+=val;
	} tint query(int l,int r=-1)const{ // [l,r]
		tint res=0; if(r==-1)r=l,l=1; --l;//(0,r]
		for(;l<r;r-=r&(-r))res+=bit[r];
		for(;r<l;l-=l&(-l))res-=bit[l];
		return res;
	}int binlift( tint vl )const{
		int x=0,nxt;//$min\{i|\sum_{j<=i}bit_j>=val\}$
		for(int exp=ilog2(n+1);exp>=0;--exp){
			nxt=x|(1<<exp);if(nxt<=n&&bit[nxt]<vl){
				vl-=bit[nxt];x=nxt; }
		} return x+1;
	} tint raw_query(int x)const{ tint res=0;for(;x;x-=x&(-x))res+=bit[x]; return res;}
};