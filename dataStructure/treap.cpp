/*
Autor: Oscar Vargas Pabon

Notar que este treap construye un max-heap en el atributo 'key'.
No tiene anadido ningun tipo de 'lazy propagation'
Las operaciones asumen indexacion en 0

Notar que en esta version, el valor del nodo (vl) funciona como
	el tamaño del subarbol que representa (trabajo sobre una lista implicita)

Its a ( max-heap _ ascendent bst )

Testeado en https://www.spoj.com/problems/GSS6/en/ ; https://codeforces.com/problemset/problem/1558/D
*/
// std::mt19937_64 rng_64( std::chrono::steady_clock::now().time_since_epoch().count() );
// asumo el anterior elemento de mi template
namespace treap{
struct Lazy{
    int x;
    Lazy(){x=0;}
    Lazy(int xx):x(xx){};
    void combine(const Lazy&o){
    	x+=o.x; // o is a 'parent'
    }
}; struct Data{
	int pref,suf,sum,mx;
	Data(){ pref=suf=mx=-1e9; sum = 0; }
	Data( int vl ) { pref=suf=sum=mx= vl; }
	Data operator +( const Data &o ) const {
		Data res; // NO necesariamente es asociativo
		res.sum = sum+o.sum;
		res.pref=std::max(pref,sum+o.pref);
		res.suf=std::max(o.suf,o.sum+suf);
		res.mx=std::max({mx,o.mx,suf+o.pref});
		return res;
	} void combine(const Lazy&o){
        mx+=o.x;
    }
}; struct Node{
	Data mv; // mv->MyValue;
	Data sv; // sv->SubtreeValue;
	Lazy lt; // lt->LazyTag
	int key, vl; //(key,vl) max-heap by key ; bst by vl
	Node *l,*r; // hijos
	Node( const Data &dt=Data() ){
		key = rng_64()%(1<<31); vl = 1;
		mv=dt;
		sv=dt;
		l=r=NULL;
	}
}; std::ostream &operator<<(std::ostream &os,const Node*rt) {
    auto aux=[&](auto rec, const Node*nd)->void{
        if(!nd)return;
        rec(rec,nd->l); os << nd->mv.sum << ' '; rec(rec,nd->r);
    }; aux(aux,rt); return os;
} inline int gvl(Node*act){return act?act->vl:0;}
inline Data gsv(Node*act){return act?act->sv:Data();}
inline void update(Node*act){if(!act)return;
	act->sv = gsv(act->l) + act->mv + gsv(act->r);
	act->vl = gvl(act->l) +     1   + gvl(act->r);
} inline void push(Node*act){if(!act)return;
    if(act->l)act->l->lz.combine(act->lz);
    if(act->r)act->r->lz.combine(act->lz);
    act->mv.combine(act->lz);
    act->sv.combine(act->lz);
    act->lz=Lazy();
} void join( Node *&t, Node *l, Node *r){
	push(l);push(r);
	if ( !l || !r ) t=(l)?l:r;
	else if ( l->key < r->key ) join(r->l,l,r->l),t=r;
	else join(l->r,l->r,r),t=l;
	update(t);
} void split( Node*t, int x, Node *&l, Node *&r ) {
	push(t);// l=t[0..x); r=t[x..n)
	int lvl=(t)?((t->l)?t->l->vl+1:1):0;
	if (!t) l=r=NULL;
	else if ( lvl <= x ) split(t->r,x-lvl,t->r,r),l=t;
	else split(t->l,x,l,t->l),r=t;
	update(t);
} void insert(Node *&t, int pos, Data vl){
	// a[0..pos]+vl+a(pos..n)
	Node *l,*r;
	split(t,pos,l,r);
	
	Node *neo= new Node(vl);
	join(t,l,neo);
	join(t,t,r);
} void update(Node *&t, int pos, Data vl){
	Node *l,*m,*r; split(t,pos,l,m); split(m,1,m,r);
	m->mv=vl;
	join(t,l,m); join(t,t,r);
} void remove(Node *&t,int pos){
	// a[0..pos)+a(pos..n)
	Node *l,*m,*r;
	split(t,pos,l,m);
	split(m,1,m,r);
	if(m) delete m;
	join(t,l,r);
} Data query(Node *&t,int ql,int qr){
	// a[ql..qr)
	Node *l,*m,*r;
	split(t,qr,m,r);
	split(m,ql,l,m);
	
	Data res=m->sv;
	join(t,m,r); join(t,l,t);
	return res;
} Node * build(const std::vector<Data> &arr){
	std::vector<Node*> ms;//MonotonicStack
	for(const Data&act:arr){
		Node *nd=new Node(act),*prv=NULL;
		while(!ms.empty()&&ms.back()->key<=nd->key){
			prv=ms.back(); ms.pop_back();
			update(prv);
		} if(!ms.empty())ms.back()->r=nd;
		nd->l=prv;
		ms.push_back(nd);
		update(nd);
	} Node *prv=NULL;
	while(!ms.empty()){
		ms.back()->r=prv;
		prv=ms.back(); ms.pop_back();
		update(prv);
	} return prv;
}
}//end namespace treap