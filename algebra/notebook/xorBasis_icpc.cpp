const int lgi=30; // this impl works entirely by noticing that
// I only have control on 'active' bits. The 'time' is assumed
struct XorB{ vector<pair<int,int>> bs;//to be increasing
	XorB(){bs=vector<pair<int,int>>(lgi,{0,-1});}
	int redu(int x,int y){//reduces element x in times
		for(int i=lgi-1;i>=0;--i)if(bs[i].second>=y)// >=y
			x=min(x,x^bs[i].first);
		return x;
	} void add(int x,int y){ // add element x in time y
		pair<int,int> act={x,y};
		for(int i=lgi-1;i>=0;--i)if((act.first>>i)&1){
			if(bs[i].second>act.second)act.first^=bs[i].first;
			else{ bs[i].first^=act.first;
				swap(bs[i],act); } }
	} // note that ord(kth(k,y),y)==k and kth(ord(x,y),y)==x
	int kth(int k,int y){
		// returns which is the kth vector (in increasing order) of
		int msk=0;for(int i=0;i<lgi;++i)if(bs[i].second>=y) {
			msk|=(k&1)<<i; k>>=1;// the span of bs(>=y)
		} int rs=0;for(int i=lgi-1;i>=0;--i)if(bs[i].second>=y){
			if(((rs>>i)&1)^((msk>>i)&1))rs^=bs[i].first;
		} return rs;
	} int ord(int x,int y){
		// returns which kth does x have on the span of bs(>=y)
		int k=0,e=0;for(int i=0;i<lgi;++i)if(bs[i].second>=y){
			k|=((x>>i)&1)<<e;
			++e;
		} return k;
	} int fred(int y){int rs=0;// degrees of freedom
		for(int i=0;i<lgi;++i)if(bs[i].second>=y)++rs;
		return rs;
	} int mx(int t){ int rs=0; // max element in span
		for(int i=lgi-1;i>=0;--i)if(tm[i]>=t)
			rs=max(rs,rs^bs[i]);
		return rs; } };