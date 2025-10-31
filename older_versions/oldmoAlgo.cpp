/*
Autor: Oscar Vargas Pabon

Funciona en tiempo

O((n\sqrt{n}+q\sqrt{n})*(T(REM)+T(ADD)))
Donde T(REM) y T(ADD) son los tiempos de remover y añadir un indice
	del global que maneja MO.

Probado con https://www.spoj.com/problems/DQUERY/

Encontre una manera de generalizarlo y no hacer estricto los 'bloques'
*/

int sqr;
class Query{
public:
    int ind, l, r;
    Query()=default;
    bool operator < ( const Query &o ) const{
		// los pares en orden decreciente
		// los impares en orden creciente
		return ((l/sqr)&1) ? r>o.r : r<o.r;
	}
};


vector<int> mo_algo( int n, const vector<Query> &arr ) {
	int q = arr.size();
	
	sqr=1; // encuentro la raiz cuadrada en tiempo \sqrt{n}
    while ( sqr*sqr < n ) ++sqr;
	
	
	
	vector<vector<Query>> q_box(sqr+1); // construyo las cajas
	for ( int i = 0 ; i < q ; ++i ) q_box[arr[i].l/sqr].push_back(arr[i]);
	for ( int i = 0 ; i <= sqr ; ++i ) sort(q_box[i].begin(),q_box[i].end());
	
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
	for ( int i = 0 ; i <= sqr ; ++i ){
        for ( const Query &act : q_box[i] ){
            while ( r > act.r ) rem(r--);
            while ( r < act.r ) add(++r);
            while ( l < act.l ) rem(l++);
            while ( l > act.l ) add(--l);
			
            res[act.ind] = cnt; // esto depende del ejercicio
        }
    }
	return res;
}
