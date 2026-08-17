/*
Autor: Oscar Vargas Pabon

pseudo-tested in https://codeforces.com/gym/106463/problem/C
				 https://codeforces.com/problemset/problem/600/E
If I did it correctly, it can be persistent, it can be sparse.
	It can be anything xdxdxd
DONT take out the 'tmp' vals. Turns out there is a weird bug in GCC in which
	the <variable>=<expression> evaluates first <variable> and if evaluating
	<expression> changes the position in memory of <variable> then it breaks
*/
template<bool persistent=0,typename tint=int>
struct UnSegTree{   
    struct Lazy{
        int vl,vsm;
        Lazy(int v=0,int vs=0):vl(v),vsm(vs){};
        
        void combine( const Lazy &l ){
            // unir con otro lazy (l estaba en el padre)
            vl+=l.vl;vsm+=l.vsm;
        } pair<Lazy,Lazy> separe( int l, int r ) {
            assert(l==l&&r==r);// separar en izquierda y derecha
            return {*this,*this};
        } bool operator == (const Lazy&o)const{return vl==o.vl;}
    }; struct Data {
        lint sum;int freq;
        Data(lint s=0,int frq=0) : sum(s),freq(frq){};
        Data operator + ( const Data &o ) const {
            Data rs(0,max<int>(freq,o.freq));
            if(rs.freq==  freq)rs.sum+=  sum;
            if(rs.freq==o.freq)rs.sum+=o.sum;
            return rs;
        } void update( const Lazy &o, int l, int r ){
            assert(l==l&&r==r);// modificar el dato
            freq+=o.vl;sum+=o.vsm;
        }
    }; struct Node { array<int,2> chld; Data dat; Lazy tag; Node(){chld={0,0};dat=Data();tag=Lazy();}; };
    int t_sz;tint rt_n; vector<Node> pool; vector<int> to_erase; inline void reset()noexcept{rt_n=1;to_erase.clear();}
    UnSegTree(tint ttn, int reser=1){t_sz=ttn-1;pool.reserve(reser);pool.pb(Node());rt_n=1;}
    inline int new_node(){
        if(!to_erase.empty()){
            int nnd;pool[nnd=to_erase.back()]=Node();to_erase.pop_back();return nnd;
        } if(int(pool.size())==rt_n)pool.push_back(Node()); pool[rt_n]=Node(); return rt_n++;}

    int build( const vector<Data> &arr, int l=0,int r=-1 ) { if(r==-1)r=t_sz;
        const int nd = new_node(); if ( l>=r ) pool[nd].dat = arr[l];
        else { const int m = (l+r)/2;
            const int tc1 = build(arr, l ,m), tc2 = build(arr,m+1,r);
            pool[nd].chld[0]=tc1;pool[nd].chld[1]=tc2;
            pool[nd].dat = pool[tc1].dat + pool[tc2].dat;
        } return nd;
    } void push(int nd,tint l,tint r) { if ( pool[nd].tag==Lazy() ) return;
        pool[nd].dat.update( pool[nd].tag,l,r ); if ( l < r ) {
            const int ln=(persistent||!pool[nd].chld[0])?new_node():pool[nd].chld[0],
                      rn=(persistent||!pool[nd].chld[1])?new_node():pool[nd].chld[1];
            if(persistent&&pool[nd].chld[0])pool[ln]=pool[pool[nd].chld[0]];
            if(persistent&&pool[nd].chld[1])pool[rn]=pool[pool[nd].chld[1]];
            // not always necessary this part of separe
            pair<Lazy,Lazy> tsep=pool[nd].tag.separe(l,r);
            pool[nd].chld[0]=ln; pool[ln].tag.combine(tsep.first );
            pool[nd].chld[1]=rn; pool[rn].tag.combine(tsep.second);
        } pool[nd].tag = Lazy(); // para no sobrecontar
    } int update( int nd,tint ql,tint qr, const Lazy &upd,tint l=0,tint r=-1){ if(r==-1)r=t_sz;
        if(!nd)nd=new_node();
        push(nd,l,r);
        // Data&mdat=pool[nd].dat;if(mdat.mx<=upd.vl)return nd;  // non-standard line (segmentTreeBeats)
        int neo = nd; if ( ql<=l && r <= qr ){//&& mdat.mx>upd.vl&& mdat.smx<upd.vl ) {
            if(persistent){neo=new_node();pool[neo]=pool[nd];}// non-standard cond (SegmentTreeBeats)
            pool[neo].tag=upd; push(neo,l,r); 
        } else if (!(r<ql||qr<l)){ const tint m =(l+r)/2; if(persistent)neo=new_node();
            const int tc1=update(pool[nd].chld[0],ql,qr,upd, l ,m);
            const int tc2=update(pool[nd].chld[1],ql,qr,upd,m+1,r);
            pool[neo].chld[0]=tc1; pool[neo].chld[1]=tc2;
            pool[neo].dat=pool[tc1].dat+pool[tc2].dat;
        } return neo;
    } Data query(int nd,tint ql,tint qr,tint l=0,tint r=-1){if(r==-1)r=t_sz;
        if(!nd)nd=new_node();
        Data res; push(nd,l,r); if ( ql<=l && r <= qr ) res = pool[nd].dat;
        else if ( !(r<ql||qr<l) ) { const tint m = (l+r)/2;
            res=query(pool[nd].chld[0],ql,qr,l,m)+query(pool[nd].chld[1],ql,qr,m+1,r);
            pool[nd].dat=pool[pool[nd].chld[0]].dat+pool[pool[nd].chld[1]].dat;
        } return res;
    } int merge(int to,int from,tint l=0,tint r=-1){if(r==-1)r=t_sz;
    // this here is tested for cases in which there is no Lazy tags
        assert(!persistent);//no such support is given
        if(!to)return from;
        else if(!from)return to;
        if(l<r){ const tint m=(l+r)/2ll;
            Node&nt=pool[to],&nf=pool[from];
            if(nt.chld[0]==nt.chld[1]&&!nt.chld[0]){
                push(from,l,r);nf.tag=nt.tag;
                nt=nf; push(to,l,r);
            }else if(nf.chld[0]==nf.chld[1]&&!nf.chld[0]){
                push(to,l,r);nt.tag=nf.tag; push(to,l,r);
            }else{ push(to,l,r);push(from,l,r);
                const array<int,2> bndl={l,m+1},bndr={m,r}; for(int i=0;i<2;++i){
                    const int tmp=merge(nt.chld[i],nf.chld[i],bndl[i],bndr[i]);
                    pool[to].chld[i]=tmp;
                } nt.dat=pool[nt.chld[0]].dat+pool[nt.chld[1]].dat;
            }
        } else{ //Depends on how to aggregate the dat values
            push(to,l,r);push(from,l,r);
            pool[to].dat.freq+=pool[from].dat.freq;
            pool[to].dat.sum=max<lint>(pool[to].dat.sum,pool[from].dat.sum);
        } to_erase.push_back(from);
        return to;
    }
}; typedef UnSegTree<bool(0)>SegT; typedef SegT::Lazy Lazy; typedef SegT::Data Data;




/////other shit https://www.luogu.com.cn/problem/P3521
struct Node{ int sm;
    Node*l,*r;
    Node()=default;
};int t_sz;
inline int get_sm(Node*nd)noexcept{return nd?nd->sm:0;}
Node*update(Node*nd,int x,int dl,int l=0,int r=-1){if(r==-1)r=t_sz;
    if(!nd)nd=new Node();
    if(l>=r)nd->sm+=dl;
    else {const int m=(l+r)/2;
        if(x<=m)nd->l=update(nd->l,x,dl,l,m);
        if(x>m)nd->r=update(nd->r,x,dl,m+1,r);
        nd->sm=get_sm(nd->l)+get_sm(nd->r);
    } return nd;
} int query(Node*nd,int ql,int qr,int l=0,int r=-1){if(r==-1)r=t_sz;
    if(!nd)return 0;
    int res=0;if(ql<=l&&r<=qr)res=nd->sm;
    else if(!(r<ql||qr<l)){ const int m=(l+r)/2;
        res=query(nd->l,ql,qr,l,m)+query(nd->r,ql,qr,m+1,r);
    } return res;
} int pquery(Node*nd,int x,int l=0,int r=-1){if(r==-1)r=t_sz;
    if(!nd)return 0;
    int res;if(l>=r)res=nd->sm;
    else { const int m=(l+r)/2;
        if(x<=m)res=pquery(nd->l,x,l,m);
        else res=pquery(nd->r,x,m+1,r);
    } return res;
} Node*merge(Node*to,Node*from,int l=0,int r=-1){if(r==-1)r=t_sz;
    if(!to)return from;
    if(!from)return to;
    if(l>=r)to->sm+=from->sm;
    else{const int m=(l+r)/2;
        to->l=merge(to->l,from->l,l,m);
        to->r=merge(to->r,from->r,m+1,r);
        to->sm=get_sm(to->l)+get_sm(to->r);
    }delete from;
    return to;
}void free_tree(Node*nd){if(nd){free_tree(nd->l),free_tree(nd->r);delete nd;}}