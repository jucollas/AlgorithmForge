/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe en testing\rangeQtest.cpp
*/
template<typename tint> struct fenw{
	vector<tint>bit; int n;
	fenw(int nn=0,int vl=0){ n=nn;bit.resize(nn);}
	fenw(const vector<tint>&arr) { assert(arr[0]==0);//note arr[0] is undefined
		n=arr.size();bit.resize(nn,0);fill(bit.begin(),bit.end(),tint(0));
		int nxt; for(int i=0;i<=n;++i) {
			bit[i]+=arr[i];nxt=i+(i&(-i));
			if(nxt<=n)bit[nxt]+=bit[i]; }
	} void update( int x, tint val ) {
		for(;x<=n;x+=x&(-x))bit[x]+=val;
	} tint query(int l,int r=-1){ // [l,r]
		tint res=0; if(r==-1)r=l,l=1; --l;//(0,r]
		for(;l<r;r-=r&(-r))res+=bit[r];
		for(;r<l;l-=l&(-l))res-=bit[l];
		return res;
	}int binlift( tint vl ) {
		int x=0,nxt;//$min\{i|\sum_{j<=i}bit_j>=val\}$
		for(int exp=ilog2(n+1);exp>=0;--exp){
			nxt=x|(1<<exp);if(nxt<=n&&bit[nxt]<vl){
				vl-=bit[nxt];x=nxt; }
		} return x+1; }
};

//////////////////// previous version
int bit[template_limit], bit_size; // BIT -> Binary Indexed Tree

void build( int *arr, int n ) {
	/* Crea el arbol en tiempo O(n) sobre el arreglo 'arr' */
	int nxt, bit_size = n+1; // porque esta indexado en 1
	memset( bit, 0, sizeof(int)*bit_size );
	for ( int i = 0 ; i < bit_size ; ++i ) {
		bit[i] += arr[i]; nxt = i + (i&(-i));
		if ( nxt < bit_size ) bit[nxt] += bit[i];
	}
}
int query( int x ) {
	/* Responde la suma del prefijo de bit[1]+...+bit[x] en tiempo O(lg n)*/
	int res = 0;
	for ( ; x > 0 ; x -= x&(-x) ) res += bit[x];
	return res;
}
int query( int l, int r ) {
	/* responde la suma del rango de bit[l]+...+bit[r] en tiempo O(lg n) */
	return query( r ) - query( l-1 );
}
int query2(int l,int r){
	// la misma query de arriba pero funciona con trucos que me saque
	// del servidor de disc GF comentario de Aeren en tips-and-tricks
	int res=0;
	for(;l<r;r-=r&(-r))res+=bit[r];
	for(;r<l;l-=l&(-l))res-=bit[l];
	return res;
}
void update( int x, int val ) {
	/* Modifica el valor de arr[x] en la representacion arborea en O(lg n) */
	for ( ; x < bit_size ; x += x&(-x) ) bit[x] += val;
}

// computa el logaritmo entero de num(la cantidad de bits activos) - esta fun es aparte
//int ilog2( int num ) { return sizeof(int)*8 - __builtin_clz( num ) -1; } 
// asumo el anterior de mi template
int binlift( int val ) {
	/* Retorna el menor indice 'x' tal que query( x ) >= val en tiempo O(lg n) */
	int x = 0, nxt;
	for ( int exp = ilog2(bit_size) ; exp >= 0 ; --exp ) {
		nxt = x|(1<<exp);
		if ( nxt < bit_size && bit[nxt] < val ) {
			val -= bit[nxt];
			x = nxt;
		}
	}
	return x+1;
}

