/*
Autor: Oscar Vargas Pabon

Notar que este treap construye un max-heap en el atributo 'key'.
No tiene anadido ningun tipo de 'lazy propagation'
Las operaciones asumen indexacion en 0

Testeado en https://www.spoj.com/problems/GSS6/en/
*/

// mt19937_64 random_64( chrono::steady_clock::now().time_since_epoch().count() );
// asumo el anterior elemento de mi template

class Data{
	public:
	int pref,suf,sum,mx;
	Data(){
		pref=suf=mx=-1e9;
		sum = 0;
	}
	Data( int vl ) { pref=suf=sum=mx= vl; }
	Data operator +( const Data &o ) const {
		Data res; // NO necesariamente es asociativo
		res.sum = sum+o.sum;
		res.pref=max(pref,sum+o.pref);
		res.suf=max(o.suf,o.sum+suf);
		res.mx=max({mx,o.mx,suf+o.pref});
		return res;
	}
};

class Node{
	public:
	// mv->MyValue; sv->SubtreeValue;
	Data mv,sv;
	int key, sz;
	Node *l,*r; // hijos
	Node( const Data &dt=Data() ){
		key = random_64()%(1ll<<30); sz = 1;
		mv=sv=dt;
		l=r=NULL;
	}
	void update(){
		auto subt_dt=[&]( Node *act ){
			return (act)?act->sv:Data();
		};
		sv = subt_dt(l) + mv + subt_dt(r) ;
		auto subt_sz=[&]( Node *act ){
			return (act)?act->sz:0;
		};
		sz = subt_sz(l) + 1 + subt_sz(r);
	}
};
void join( Node *&t, Node *l, Node *r){
	if ( !l || !r ) t=(l)?l:r;
	else if ( l->key < r->key ) join(r->l,l,r->l),t=r;
	else join(l->r,l->r,r),t=l;
	t->update();
}
void split( Node* t, int x, Node *&l, Node *&r ) {
	// l=t[0..x); r=t[x..n)
	int lsz=(t)?((t->l)?t->l->sz+1:1):0;
	if (!t) l=r=NULL;
	else if ( lsz <= x ) split(t->r,x-lsz,t->r,r),l=t;
	else split(t->l,x,l,t->l),r=t;
	if(t)t->update();
}

void t_in(Node *&t, int pos, Data vl){
	// a[0..pos)+vl+a[pos..n)
	Node *l,*r;
	split(t,pos,l,r);
	
	Node *neo= new Node(vl);
	join(t,l,neo);
	join(t,t,r);
}

void t_upd(Node *&t, int pos, Data vl){
	Node *l,*m,*r; split(t,pos,l,m); split(m,1,m,r);
	m->mv=vl;
	join(t,l,m); join(t,t,r);
}

void t_out(Node *&t,int pos){
	// a[0..pos)+a(pos..n)
	Node *l,*m,*r;
	split(t,pos,l,m);
	split(m,1,m,r);
	if(m) delete m;
	join(t,l,r);
}

Data t_query(Node *&t,int ql,int qr){
	// a[ql..qr)
	Node *l,*m,*r;
	split(t,qr,m,r);
	split(m,ql,l,m);
	
	Data res=m->sv;
	join(t,m,r); join(t,l,t);
	return res;
}

Node * t_build(const vector<Data> &arr){
	vector<Node*> ms;//MonotonicStack
	for(const Data &act:arr){
		Node *nd=new Node(act),*prv=NULL;
		while(!ms.empty()&&ms.back()->key<=nd->key){
			prv=ms.back(); ms.pob();
			prv->update();
		}
		if(!ms.empty())ms.back()->r=nd;
		nd->l=prv;
		ms.pb(nd);
		nd->update();
	}
	Node *prv=NULL;
	while(!ms.empty()){
		ms.back()->r=prv;
		prv=ms.back(); ms.pob();
		prv->update();
	}
	return prv;
}

void t_debug(Node*t){
	function<void(Node*)> aux=[&](Node*nd){
		if(!nd)return;
		aux(nd->l);
		
		cout << "( " << nd->key << " ";
		cout << nd->sv.sum << ' ';
	
		aux(nd->r);
	}; aux(t); cout << endl;
}