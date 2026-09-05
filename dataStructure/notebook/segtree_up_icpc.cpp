struct Data{
    int mn,mx;
    Data(int n=inf,int x=-inf):mn(n),mx(x){};
    Data operator +(const Data&o)const{
        return Data( min<int>(mn,o.mn),max<int>(mx,o.mx) ); }
// remember I can modify it to receive range-update and point-query
}; struct SegT{// in a fairly straightforward manner
	vector<Data> tree;int n;
	SegT()=default;
	SegT(int nn):n(nn){ tree.resize(2*nn); }
	SegT(const vector<Data>&arr){
		n=arr.size();tree.resize(2*n);//O(n)
		for(int i=0;i<n;++i)tree[i+n]=arr[i];
		for(int i=n-1;i;--i)tree[i]=tree[i*2]+tree[i*2+1];
	} void update( int x,const Data&val ) {
	/* Modifica el valor en arr[x] segun la representacion del arbol en O(lg n) */
		tree[x+=n]=val;
		for(x/=2;x;x/=2)tree[x]=tree[x*2]+tree[x*2+1];
	} Data query( int l, int r )const{
	/* Responde a la query arr[l] + ... + arr[r] en tiempo O(lg n) */
		Data res;
		for (l+=n,r+=n;l<=r;l/=2,r/=2) {
			if(  l&1 )res=res+tree[l++];
			if(!(r&1))res=res+tree[r--];
		} return res; } };