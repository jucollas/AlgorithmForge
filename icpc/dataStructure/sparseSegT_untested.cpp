
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

class Node {
	public:
	Node *l,*r;
	Data dat; Lazy tag;
	Node(){ l=r=NULL; }
};
int t_sz=0; Node*rt=NULL;
void push( Node*nd, int l, int r ) {
/* Empuja el valor del tag a los nodos mas profundos */
	nd->dat.update( nd->tag,l,r );
	if ( l < r ) {
		if ( !nd->l ) nd->l = new Node();
		if ( !nd->r ) nd->r = new Node();
		
		pair<Lazy,Lazy> tsep=nd->tag.separe(l,r);
		nd->l->tag.combine(tsep.first);
		nd->r->tag.combine(tsep.second);
	}
	nd->tag = Lazy(); // para no sobrecontar
}
Data query(int ql, int qr, Node*nd=rt,int l=0,int r=t_sz){
	Data res; push(nd,l,r);
	if( ql<=l&&r<=qr ) res=nd->dat;
	else if (!(r<ql||qr<l)){
		int m=(l+r)/2;
		res = query(ql,qr,nd->l,l,m)+query(ql,qr,nd->r,m+1,r);
	}
	return res;
}
void update(int ql, int qr, const Lazy &upd,Node*nd=rt,int l=0, int r=t_sz){
	push(nd,l,r);
	if( ql<=l&&r<=qr ) {
		nd->tag = upd; // porque ya hize el push antes
		push(nd,l,r);
	} else if (!(r<ql||qr<l)){
		int m=(l+r)/2;
		update(ql,qr,upd,nd->l,l,m);
		update(ql,qr,upd,nd->r,m+1,r);
		nd->dat=nd->l->dat+nd->r->dat;
	}
	return nd;
}


void build(int n){
	t_sz=n-1;
	rt=new Node();
	// rt->l=rt->r=rt; // no entiendo el truquito del santo
}