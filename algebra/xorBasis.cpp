/*
Autor: Oscar Vargas Pabon

Remember that the amount of ways of generating x, its 2**(|S|-|bs|) if x is generable
*/

// ultra slim impl
// untested
int red(const vector<int>&bs,int x){for(int ac:bs)x=min(x,x^ac);return x;}
bool add(vector<int>&bs,int x){x=red(x);if(x)bs.pb(x);return x;}
int mx(const vector<int>&bs){int x=0;for(int ac:bs)x=max(x,x^ac);return x;}

// esta idea la tuve en la superregional
// no esta testeado
class XorR{
	public:
	vector<int> bs,us,ori;
	XorR()=default;
	int red(int x){for(int ac:bs)x=min(x,x^ac);}
	int mx(){int x=0;for(int ac:bs)x=max(x,x^ac);}
	
	vector<int> raw_rec(int msk){
		vector<int> res;
		rep(i,0,n)if(msk&(1<<i))res.pb(ori[i]);
		return res;
	}
	vector<int> rec(int x){
		// reconstruye los elementos que generan x
		int msk=0;
		rep(i,0,bs.size())if(x>(x^bs[i])){
			x^=bs[i]; msk^=us[i];
		}
		return raw_rec(msk);
	}
	vector<int> mx_rec(){
		// reconstruye los elementos que generan el maximo
		int msk=0,x=0;
		rep(i,0,bs.size())if(x<(x^bs[i])){
			x^=bs[i]; msk^=us[i];
		}
		return raw_rec(msk);
	}
	bool add(int x){
		int s=1<<bs.size(),ox=x;
		rep(i,0,bs.size())if(x>(x^bs[i])){
			x^=bs[i]; s^=us[i];
		}
		if(x)bs.pb(x),us.pb(s),ori.pb(ox);
		return x;
	}
};

// range static Xor basis
// tested in https://codeforces.com/contest/1100/problem/F
const int lgi=20;
class XorP{
	public:
	array<int,lgi> bs,tm;
	XorP(){
		rep(i,0,lgi)bs[i]=0;
		rep(i,0,lgi)tm[i]=-1;
	}
	
	int red(int x,int t){
		rep(i,lgi-1,-1)if(tm[i]>=t)x=min(x,x^bs[i]);
		return x;
	}
	void add(int x, int t){
		rep(i,lgi-1,-1) if( (x>>i)&1 ) {
			if(tm[i] < t ){
				swap(tm[i],t);swap(bs[i],x);
			}
			
			x^=bs[i];
		}
	}
	int mx(int t){
		int rs=0;
		rep(i,lgi-1,-1)if(tm[i]>=t)rs=max(rs,rs^bs[i]);
		return rs;
	}
};