/*
Autor: Oscar Vargas Pabon

Works in O(n lgn)
Think on impl with vector or smthng

Tested in https://codeforces.com/problemset/problem/1017/E
*/

vector<Pt> convex_hull(Pt *arr,int n){
	sort(arr,arr+n);
	vector<Pt> up,dwn;
	rep(i,0,n){
		Pt &ac=arr[i]; int pz=up.size(),dz=dwn.size();
		while(pz>1&& (ac-up[pz-1])%(up[pz-2]-up[pz-1])>=0 ){
			--pz; up.pob();
		}
		if(pz>=1&&up.back().x==ac.x&&up.back().y<ac.y)up.pob();
		if(up.empty()||up.back().x<ac.x||up.back().y<ac.y)up.pb(ac);
		
		while(dz>1&& (ac-dwn[dz-1])%(dwn[dz-2]-dwn[dz-1])<=0 ){
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