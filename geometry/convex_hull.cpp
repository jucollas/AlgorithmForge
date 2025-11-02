/*
Autor: Oscar Vargas Pabon

Works in O(n lgn)
Think on impl with vector or smthng
There are stuff / versions that are incremental

Tested in https://codeforces.com/problemset/problem/1017/E
          https://open.kattis.com/problems/convexhull

I assume from my template:
#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
#define pb push_back
#define pob pop_back
*/

vector<Pt> convex_hull(Pt *arr,int n){
	sort(arr,arr+n);
	vector<Pt> up,dwn;
	auto crss=[](const Pt&a,const Pt &b,const Pt &bs){return (a-bs)%(b-bs);};
	rep(i,0,n){
		Pt &ac=arr[i]; int pz=up.size(),dz=dwn.size();
		while(pz>1&& crss(ac,up[pz-2],up[pz-1])>=0 ){
			--pz; up.pob();
		}
		if(pz>=1&&up.back().x==ac.x&&up.back().y<ac.y)up.pob();
		if(up.empty()||up.back().x<ac.x||up.back().y<ac.y)up.pb(ac);
		
		while(dz>1&& crss(ac,dwn[dz-2],dwn[dz-1])<=0 ){
			--dz; dwn.pob();
		}
		if(dz>=1&&dwn.back().x==ac.x&&dwn.back().y>ac.y)dwn.pob();
		if(dwn.empty()||dwn.back().x<ac.x||dwn.back().y>ac.y)dwn.pb(ac);
	}
	vector<Pt> res=dwn;
	rep(i,up.size()-1,-1){
		if( !(res.back()==up[i]) ) res.pb(up[i]);
	}
	// idebug(up);idebug(dwn);
	if(int(res.size())>1&&res.front()==res.back())res.pob();
	return res;
}


vector<int> hullInd(const vector<Pt>& v) {
	//adapted from https://github.com/bqi343/cp-notebook
	
	int ind = int(min_element(all(v))-begin(v));
	vector<int> cand, C={ind};
	rep(i,0,v.size()) if (v[i] != v[ind]) cand.pb(i);
	sort(all(cand),[&](int a, int b) { 
		// sort by angle, tiebreak by distance
		Pt x = v[a]-v[ind], y = v[b]-v[ind]; Tpt t = x%y;
		return t != 0 ? t > 0 : x*x < y*y;
	}); 
	auto crss=[](const Pt&a,const Pt &b,const Pt &bs){return (a-bs)%(b-bs);};
	for (int c:cand){
		while (int(C.size()) > 1 && crss(v[C.back()],v[c],v[C[C.size()-2]]) <= 0)
			C.pob();
		C.pb(c); }
	return C;
}
vector<Pt> hull(const vector<Pt>& v) {
	vector<int> w = hullInd(v); vector<Pt> res;
	for(int t:w) res.pb(v[t]);
	return res; }