/*
Autor: Oscar Vargas Pabon

pseudo-tested in https://codeforces.com/gym/106463/problem/C
If I did it correctly, it can be persistent, it can be sparse.
	It can be anything xdxdxd
DONT take out the 'tmp' vals. Turns out there is a weird bug in GCC in which
	the <variable>=<expression> evaluates first <variable> and if evaluating
	<expression> changes the position in memory of <variable> then it breaks
*/
template<bool persistent=0>
struct UnSegTree{	
	struct Lazy{
		int vl;
		Lazy(int v=1e9):vl(v){};
		
		void combine( const Lazy &l ){
			// unir con otro lazy (l estaba en el padre)
			vl=min<int>(l.vl,vl);
		} pair<Lazy,Lazy> separe( int l, int r ) {
			assert(l==l&&r==r);// separar en izquierda y derecha
			return {*this,*this};
		} bool operator == (const Lazy&o)const{return vl==o.vl;}
	}; struct Data {
		lint sum;int mx,smx,cnt;
		Data(lint s=0,int mmx=-1,int ssmx=-1,int ccnt=0) : sum(s),mx(mmx),smx(ssmx),cnt(ccnt) {};
		Data operator + ( const Data &o ) const {
		    Data rs;rs.sum=sum+o.sum; rs.mx=max<int>(mx,o.mx);
		    rs.cnt=0;if(mx>=o.mx)rs.cnt+=cnt;
		    if(mx<=o.mx)rs.cnt+=o.cnt;
		    rs.smx=max<int>(smx,o.smx);if(mx!=o.mx)rs.smx=max<int>(rs.smx,min<int>(mx,o.mx));
		    return rs;
		} void update( const Lazy &o, int l, int r ){
			assert(l==l&&r==r);// modificar el dato
			if( mx<=o.vl || smx>=o.vl)return;
			sum-=(mx-o.vl)*1ll*cnt; mx=o.vl;
		}
	}; struct Node { array<int,2> chld; Data dat; Lazy tag; Node(){chld={-1,-1};}; };
	int t_sz, rt_n; vector<Node> pool; inline void reset()noexcept{rt_n=0;}
	UnSegTree(int ttn, int reser=0){t_sz=ttn-1;pool.reserve(reser);rt_n=0;}
	inline int new_node(){ if(int(pool.size())==rt_n)pool.push_back(Node()); pool[rt_n]=Node(); return rt_n++;}

	int build( const vector<Data> &arr, int l=0,int r=-1 ) { if(r==-1)r=t_sz;
		const int nd = new_node(); if ( l>=r ) pool[nd].dat = arr[l];
		else { const int m = (l+r)/2;
			const int tc1 = build(arr, l ,m), tc2 = build(arr,m+1,r);
			pool[nd].chld[0]=tc1;pool[nd].chld[1]=tc2;
			pool[nd].dat = pool[tc1].dat + pool[tc2].dat;
		} return nd;
	} void push( int nd, int l, int r ) { if ( pool[nd].tag==Lazy() ) return;
		pool[nd].dat.update( pool[nd].tag,l,r ); if ( l < r ) {
			const int ln=(persistent||pool[nd].chld[0]==-1)?new_node():pool[nd].chld[0],
				      rn=(persistent||pool[nd].chld[1]==-1)?new_node():pool[nd].chld[1];
			if(persistent&&pool[nd].chld[0]!=-1)pool[ln]=pool[pool[nd].chld[0]];
			if(persistent&&pool[nd].chld[1]!=-1)pool[rn]=pool[pool[nd].chld[1]];
			// not always necessary this part of separe
			pair<Lazy,Lazy> tsep=pool[nd].tag.separe(l,r);
			pool[nd].chld[0]=ln; pool[ln].tag.combine(tsep.first );
			pool[nd].chld[1]=rn; pool[rn].tag.combine(tsep.second);
		} pool[nd].tag = Lazy(); // para no sobrecontar
	} int update( int nd, int ql,int qr, const Lazy &upd, int l=0, int r=-1 ){ if(r==-1)r=t_sz;
		if(nd==-1)nd=new_node();
		push(nd,l,r);
		Data&mdat=pool[nd].dat;if(mdat.mx<=upd.vl)return nd;  // non-standard line (segmentTreeBeats)
		int neo = nd; if ( ql<=l && r <= qr && mdat.mx>upd.vl&& mdat.smx<upd.vl ) {
		    if(persistent){neo=new_node();pool[neo]=pool[nd];}// non-standard cond (SegmentTreeBeats)
			pool[neo].tag=upd; push(neo,l,r); 
		} else if (!(r<ql||qr<l)){ const int m =(l+r)/2; if(persistent)neo=new_node();
			const int tc1=update(pool[nd].chld[0],ql,qr,upd, l ,m);
			const int tc2=update(pool[nd].chld[1],ql,qr,upd,m+1,r);
			pool[neo].chld[0]=tc1; pool[neo].chld[1]=tc2;
			pool[neo].dat=pool[tc1].dat+pool[tc2].dat;
		} return neo;
	} Data query( int root, int ql, int qr ) {
		function<pair<int,Data>(int,int,int)>raw_query=[&](int nd,int l,int r)->pair<int,Data>{
			if(nd==-1)nd=new_node();
			Data res; push(nd,l,r); if ( ql<=l && r <= qr ) res = pool[nd].dat;
			else if ( !(r<ql||qr<l) ) { const int m = (l+r)/2;
				pair<int,Data> lc=raw_query(pool[nd].chld[0],l,m),rc=raw_query(pool[nd].chld[1],m+1,r);
				res=lc.second+rc.second;
				pool[nd].chld[0]=lc.first; pool[nd].chld[1]=rc.first;
			} return {nd,res};
		}; return raw_query( root,0,t_sz ).second;		
	} int merge(int to,int from,int l=0,int r=-1){if(r==-1)r=t_sz;
		// NOTE: this operation is currently untested
		assert(!persistent);//no such support is given
		if(to==-1)return from;
		else if(from==-1)return to;
		push(to,l,r);push(from,l,r);
		to=update(to,pool[from].dat,l,r); const lint m=(l+r)/2ll;
		push(to,l,r);
		const array<lint,2> bndl={l,m+1},bndr={m,r}; for(int i=0;i<2;++i){
			const int tmp=merge(pool[to].chld[i],pool[from].chld[i],bndl[i],bndr[i]);
			pool[to].chld[i]=tmp;
		} return to;
	}
}; typedef UnSegTree<bool(1)>PSTree; typedef PSTree::Lazy Lazy; typedef PSTree::Data Data;