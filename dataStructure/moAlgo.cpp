/*
Autor: Oscar Vargas Pabon

Recordar que existe la idea de balanceo (usar sqrt-decomp o otra de modo que 
	la update se haga en O(1), aunque la query pueda empeorar, pues el algoritmo hace
	mas updates que queries)


Asumo que T(REM) y T(ADD) son los tiempos de remover y anadir un indice
	del global que maneja MO.

Idea 1: usar bloque de tamano $\sqrt(n)$. Genera tiempo $O((n+q)\sqrt(n))$
Idea 2: usar bloque de tamano $\frac{n}{\sqrt(q)}$. Genera tiempo $O(n\sqrt(q))$
	source: https://codeforces.com/blog/entry/61203?#comment-451304
Idea 3: Usar orden de hilbert. Genera tiempo $O(n\sqrt(q))$
	source: https://codeforces.com/blog/entry/61203
	impl:https://codeforces.com/blog/entry/61203?#comment-1064868

Nota: La idea 1 (o tal vez 2 tambien) puede no funcionar muy bien cuando se requiere hacer un balanceo
Probado con https://www.spoj.com/problems/DQUERY/
Otro problema interesante https://codeforces.com/problemset/problem/2006/D
*/
uint64_t hilbertorder(uint64_t x, uint64_t y) {
	// https://codeforces.com/blog/entry/61203?#comment-1064868
    const uint64_t logn = __lg(max(x, y) * 2 + 1) | 1;
    const uint64_t maxn = (1ull << logn) - 1;
    uint64_t res = 0;
    for (uint64_t s = 1ull << (logn - 1); s; s >>= 1) {
        bool rx = x & s, ry = y & s;
        res = (res << 2) | (rx ? ry ? 2 : 1 : ry ? 3 : 0);
        if (!rx) {
            if (ry) x ^= maxn, y ^= maxn;
            swap(x, y);
        }
    }
    return res;
}

class Query{
public:
    int ind, l, r, ord;
    Query()=default;
    bool operator < ( const Query &o ) const{ return ord<o.ord; }
};

const int method=0;
vector<int> mo_algo( int n, vector<Query> &query ) {
	int q = query.size();
	
	// idea 1 -> $O(n\sqrt{q})$
	if(method==0)for(Query &q:query)q.ord=hilbertorder(q.l,q.r);
	else{
	  int block_size;
	  if(method==1){
		// idea 2 -> $O(n\sqrt{q})$
	    int sqrt=1; while(sqrt*sqrt<q)++sqrt;
	    block_size=n/sqrt;
	  } else {
	    // idea 3 ->$ O((n+q)\sqrt{n})$
	    block_size=1;
	    while ( block_size*block_size < n ) ++block_size;
	  }
	  block_size=max(block_size,1);
	  for(Query &q:query){
	    int bl_ind=q.l/block_size;
	    q.ord=(bl_ind&1)?n-q.r:q.r;
	    q.ord+=bl_ind*n;
	  }
	}
	
	sort(query.begin(),query.end());
	
	/// Llenar esta parte con lo necesario del problema
	vector<int> res(q);
	vector<int> all(n,0); 
	int cnt=0;
	function<void(int)> add = [&]( int x ){
		if ( !all[a[x]] ) ++cnt;
		++all[a[x]];
	};
	function<void(int)> rem = [&]( int x ) {
		--all[a[x]];
		if ( !all[a[x]] ) --cnt;
	};
	//// termina la parte del llenado
	
	int l = 0, r = -1;
    for ( const Query &act : query ){
        while ( r > act.r ) rem(r--);
        while ( r < act.r ) add(++r);
        while ( l < act.l ) rem(l++);
		while ( l > act.l ) add(--l);
		
		res[act.ind] = cnt; // esto depende del ejercicio
	}
	return res;
}
