#include<bits/stdc++.h>
using namespace std;

int main(){
	pair<int,int> ar(1,2);
	// auto [x,y]=ar;
	// cout << x << endl; cout << y << endl;
	vector<pair<int,int>> vp={{1,2},{3,4},{5,6}};
	for(auto [x,y]:vp)cout << x << " _ " << y << endl;
	
	return 0;
}