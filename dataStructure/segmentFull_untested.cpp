/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC


l -> left tree bound, r -> right tree bound
ql -> left op bound, qr -> right op bound
*/
class Lazy{
public:
	int vl;
	Lazy(int v=0):vl(v){};
	
	void combine( const Lazy &l ){
		// unir con otro lazy (l estaba en el padre)
		vl+=l.vl+vl;
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
		sum+=(r-l+1)*vl;
	}
};

class SegT{
public:
	vector<Data> tree;
	vector<Lazy> tag;
	int t_size;
	void calc_cons( int ind, int l, int r, int &m, int &rind ) {
		m = (l+r)/2; rind = ind+2*(m-l+1);
	}
void raw_build ( const vector<Data> &arr,int ind=0,int l=0,int r=t_size ) {
    /* funcion recursiva auxiliar que construye el segTree */
	
	tag[ind] = Lazy(); // siempre en neutro
	if ( l == r ) tree[ind] = arr[l]; // caso base
	else {
		int m,rind; calc_cons(ind,l,r,m,rind);
		raw_build( arr, ind+1, l, m ); // hijo izquierdo
		raw_build( arr, rind, m+1, r ); // hijo derecho
		
		tree[ind] = tree[ind+1] + tree[rind];
	}
}
void push( int ind, int l, int r ) {
/* Empuja el valor del tag a los nodos mas profundos */
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
		arr_size=arr.size();
		tag.resize(arr_size*2); tree.resize(arr.size*2);
		--arr_size;
		raw_build( arr );
	}
	

void update( int ql, int qr, Data vl, int ind=0, int l=0, int r=arr_size ) {
/* actualizo los valores en arr[lx],...,arr[rx] por arr[lx]+=val,...,arr[rx]+=val
	en el arbol en tiempo O(lg n) */
	push( ind, l, r );
	if ( ql <= l && r <= qr ) {
		tag[ind] = tag[ind] + vl;
		push( ind, l, r );
	} else if ( !(r<ql||qr<l) ) {
		int m,rind; calc_cons(ind,l,r,m,rind);
		
		update( ql, qr, val, ind+1, l, m ); // hijo izquierdo
		update( ql, qr, val, rind, m+1, r ); // hijo derecho
		tree[ind] = tree[ind+1] + tree[rind]; // actualizo el arbol
	}
}

Data query( int ql, int qr, int ind=0, int l=0, int r=arr_size ) {
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