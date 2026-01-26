/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC

tested in https://codeforces.com/problemset/problem/2117/H

l -> left tree bound, r -> right tree bound
ql -> left op bound, qr -> right op bound

For historic sums -> queries:
	1 -> A[i]+=x for l<=i<=r
	2 -> B[i]+=A[i] for 0<=i<n
	3 -> \sum_{i=l}^r B[i]
Do stuff with timestamps. Both lazy and data have 'tot,acum,time'
	normalizing to the actual time is tot'= acum*(ac_time-time)+tot
	Lazy keeps the value of a single element in acum while Data keeps the whole range
	so update is tot'=acum*(ac_time-time)+tot +(o.acum*(ac_time-o.time)+o.tot)*(r-l+1)
Idea applied in subtask of https://www.luogu.com.cn/problem/P11659
	
*/
class Lazy{
public:
	int vl;
	Lazy(int v=0):vl(v){};
	
	void combine( const Lazy &l ){
		// unir con otro lazy (l estaba en el padre)
		vl+=l.vl;
	}
	pair<Lazy,Lazy> separe( int l, int r ) {
		// separar en izquierda y derecha
		return {*this,*this};
	}
};
class Data {
	public:
	int sum;
	Data(int s=0) : sum(s) {};
	Data operator + ( const Data &o ) const {
		return Data(sum+o.sum);
	}
	void update( const Lazy &o, int l, int r ){
		// modificar el dato
		sum+=(r-l+1)*o.vl;
	}
};

class SegT{
public:
	vector<Data> tree;
	vector<Lazy> tag;
	int t_size;
	inline void calc_cons( int ind, int l, int r, int &m, int &rind ) {
		m = (l+r)/2; rind = ind+2*(m-l+1);
	}
void raw_build ( const vector<Data> &arr,int ind=0,int l=0,int r=-1 ) {
	if(r==-1)r=t_size;
    /* funcion recursiva auxiliar que construye el segTree */
	
	// tag[ind] = Lazy(); // siempre en neutro
	if ( l == r ) tree[ind] = arr[l]; // caso base
	else {
		int m,rind; calc_cons(ind,l,r,m,rind);
		raw_build( arr, ind+1, l, m ); // hijo izquierdo
		raw_build( arr, rind, m+1, r ); // hijo derecho
		
		tree[ind] = tree[ind+1] + tree[rind];
	}
}
inline void push( int ind, int l, int r ) {
/* Empuja el vlor del tag a los nodos mas profundos */
	tree[ind].update( tag[ind],l,r );
	if ( l < r ) { // no hay edge-case cuando lt==rt
		int m, rind; calc_cons(ind,l,r,m,rind);
		pair<Lazy,Lazy> sd=tag[ind].separe(l,r);
		tag[ind+1].combine(sd.first);
		tag[rind].combine(sd.second);
	}
	tag[ind] = Lazy(); // para no sobrecontar
}
	SegT()=default;
	SegT( const vector<Data>&arr ){
		t_size=arr.size();
		tag.resize(t_size*2); tree.resize(t_size*2);
		--t_size;
		raw_build( arr );
	}
	

void update( int ql, int qr, Lazy vl, int ind=0, int l=0, int r=-1 ) {
	if(r==-1)r=t_size;
/* actualizo los vlores en arr[lx],...,arr[rx] por arr[lx]+=vl,...,arr[rx]+=vl
	en el arbol en tiempo O(lg n) */
	push( ind, l, r );
	if ( ql <= l && r <= qr ) {
		tag[ind].combine(vl);
		push( ind, l, r );
	} else if ( !(r<ql||qr<l) ) {
		int m,rind; calc_cons(ind,l,r,m,rind);
		
		update( ql, qr, vl, ind+1, l, m ); // hijo izquierdo
		update( ql, qr, vl, rind, m+1, r ); // hijo derecho
		tree[ind] = tree[ind+1] + tree[rind]; // actualizo el arbol
	}
}

Data query( int ql, int qr, int ind=0, int l=0, int r=-1 ) {
	if(r==-1)r=t_size;
/* Hago una query de arr[lx]+...+arr[rx] en tiempo O(lg n) */
	Data res;
	push( ind, l, r );
	if ( ql <= l && r <= qr ) res = tree[ind];
	else if ( !(r<ql||qr<l) ) {
		int m,rind; calc_cons(ind,l,r,m,rind);
		res = query(ql,qr,ind+1,l,m)+query(ql,qr,rind,m+1,r);
	}
	return res;
}
};