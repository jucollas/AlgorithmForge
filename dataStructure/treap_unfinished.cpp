/*
xd 
*/

mt19937 random( 121651631 );

class Data{
	public:
	int val,pref,suf,sum,mx;
	Data(){
		pref=suf=mx=-1e9;
		sum = 0; val=-1;
	}
	Data( int vl ) { val=pref=suf=sum=mx= vl; }
	bool operator +( const Data &o ) const {
		Data res; // NO necesariamente es asociativo
		res.sum = sum+o.sum;
		res.pref=max(pref,sum+o.pref);
		res.suf=max(o.suf,o.sum+suf);
		res.mx=max(max(mx,o.mx),suf+o.pref);
		return res;
	}
};

class Node{
	public:
	Data dt; // datos
	int key, sz;
	Node *l,*r; // hijos
	Node( const Data &dt=Data() ){
		key = random(); sz = 1;
		this->dt = dt;
		l=r=NULL;
	}
	void update(){
		int vl = dt.val;
		auto subt =[&]( Node *act ){
			return (act)?act->dt:Data();
		};
		dt = subt(l) + Data(vl) + subt(r);
		dt.val = vl;
	}
};

void join( Node &*t, Node *l, Node *r){
	if ( !l || !r ) t=(l)?l:r;
	else if ( l->key < r->key ) join(r->l,l,r->l);
	else join(l->r,l->r,r);
	t->update();
}
void split( Node* t, int x, Node &*l, Node &*r ) {
	if (!t) l=r=NULL;
	else if ( t->sz <= x ) split(t->r,x,t->r,r),l=t;
	else split(t->l,x-t->sz,l,t->l),r=t;
	t->update();
}

void q_del( Node &*t, int x ) {
	Node *l,*m,*r,*tmp;
	split(t,x,tmp,r);
	split(tmp,x-1,l,m);
	delete m;
	join(t,l,r);
}
void q_insert( Node &*t, int x, int y ) {
	Node *l,*r,*m=new Node(y);
	split(t,x-1,l,r);
	join( t, l, m );
	join( t, t, r );
}
int q_quer( Node &*t, int x, int y ) {
	Node *l, *r, *m;
	split(t,x-1,l,m);
	split(m,y,m,r);
	int res = m->dt.mx;
	join(t,l,m);
	join(t,t,r);
}
void q_upd( Node &*t, int x, int y ) {
	Node *l, *r, *m;
	split(t,x,m,r);
	split(m,x-1,l,m);
	m->dt = Data(y);
	join(t,l,m);
	join(t,t,r);
}


