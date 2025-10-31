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


// range static XorBasis finding the kth generable element and the 
//  order of a given generable element
// used in https://codeforces.com/contest/2143/problem/F
const int lgi=30;
class XorB{
	public:
	
	vector<pair<int,int>> bs;
	XorB(){bs=vector<pair<int,int>>(lgi,{0,-1});}
	int redu(int x,int y){rep(i,lgi-1,-1)if(bs[i].second>=y)x=min(x,x^bs[i].first); return x;}
	void add(int x,int y){
		pair<int,int> act={x,y};
		rep(i,lgi-1,-1)if((act.first>>i)&1){
			if(bs[i].second>act.second)act.first^=bs[i].first;
			else{
				bs[i].first^=act.first;
				swap(bs[i],act);
			}
		}
	}
	// note that ord(kth(k,y),y)==k and kth(ord(x,y),y)==x
	int kth(int k,int y){
		// returns which is the kth vector (in increasing order) of the span of bs(>=y)
		int msk=0;rep(i,0,lgi)if(bs[i].second>=y) {
			msk|=(k&1)<<i; k>>=1;
		}
		
		int rs=0;rep(i,lgi-1,-1)if(bs[i].second>=y){
			if(((rs>>i)&1)^((msk>>i)&1))rs^=bs[i].first;
		}
		return rs;
	}
	int ord(int x,int y){
		// returns which kth does x have on the span of bs(>=y)
		int k=0,e=0;rep(i,0,lgi)if(bs[i].second>=y){
			k|=((x>>i)&1)<<e;
			++e;
		}
		return k;
	}
	
	// degrees of freedom
	int fred(int y){int rs=0;rep(i,0,lgi)if(bs[i].second>=y)++rs;return rs;}
	
};

/* Xor basis modulo m */
/* taken from errogorn's blog, maroonrk's comment https://codeforces.com/blog/entry/98376 
	I simply adapted it to my style. It seems interesting to find the kth generable elements. However,
		it may scale up quickly for certain values of mod.
*/
lint gcd(lint a, lint b, lint& x, lint& y) {
	// a*x+b*y==gcd(a,b);
    if (b == 0ll) {
        x = 1;
        y = 0;
        return a;
    }
    lint x1, y1;
    lint d = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}


const lint mod=360; const int v_sz=10;
class MBasis{
	public:
	vector<vector<lint>> bas;
	MBasis(){ bas=vector<vector<lint>>(v_sz,vector<int>(v_sz,0)); }
	void v_sum( vector<lint>&va,const vector<lint> &vb,int k){ rep(i,0,v_sz)va[i]+=vb[i]*k; }
	void add( vector<lint> &v ) {
		rep(i,0,v_sz)if(v[i]){
			int x,y;int g=gcd(v[i],bas[i][i],x,y);
			int z=(bas[i][i])?bas[i][i]/g : mod/g,w=- v[i]/g;
			rep(j,0,v_sz)tie(bas[i][j],v[j])=make_pair((v[j]*x+bas[i][j]*y)%mod,(z*v[j]+w*bas[i][j])%mod);
		}
	}
	bool red( vector<lint> &v ){
		rep(i,0,v_sz)if(v[i]){
			if(!bas[i][i]||v[i]%bas[i][i])return 0;
			v_sum(v,bas[i],- v[i]/bas[i][i]);
		}
		return 1;
	}
	
	vector<lint> mx_span(){
		vector<lint> rs(v_sz,0);
		rep(i,0,v_sz)if(bas[i][i]) v_sum( rs, bas[i], (mod-1-v[i])/bas[i][i] );
		return rs;
	}
	
	lint sz_span(){
		lint res=1;rep(i,0,v_sz)if(bas[i][i]){
			res*=mod/bas[i][i];
		}
		return res;
	}
	
};