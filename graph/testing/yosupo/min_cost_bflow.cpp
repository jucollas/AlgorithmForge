/*
 ________
|    ___ |
|  ,',.(`|
| :  `'  |
| :) _  (|
|  `:_)_,|
|________|

Autor: Oscar Vargas Pabon
Fecha: Oct 14 2026
https://judge.yosupo.jp/problem/min_cost_b_flow

NOTE: unfinished. This only works for min-cost max-flow.
*/
//#pragma GCC optimize("O3")
//#define NDEBUG
#include <bits/stdc++.h>
#include <cassert>

typedef long long lint;
// __uint128_t, __int128_t, uint64_t, int64_t, uint32_t,int32_t
using namespace std;
#ifdef OSVARP
    #include<sys/resource.h>
#else
    #define cerr for(;false;) cerr
#endif
template <typename t1,typename t2> std::ostream &operator<<(ostream &os, const std::pair<t1,t2>&pr){return os<<"("<<pr.first<<";"<<pr.second<<")";};
template <typename T> std::ostream &operator<<(ostream &os, const std::vector<T>&vc){for(auto ac:vc)os<<ac<<' ';return os;};
template <typename T> std::ostream &operator<<(ostream &os, const std::set<   T>&vc){for(auto ac:vc)os<<ac<<' ';return os;};
template <typename t1,typename t2> std::ostream &operator<<(ostream &os, const std::map<t1,t2>&vc){for(auto ac:vc)os<<ac<<' ';return os;};
#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;assert(it==it);}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define adebug(ar,n) {cerr<<'['<<#ar<<']';for(int my_imp_ind=0;my_imp_ind<n;++my_imp_ind)cerr << ' ' << ar[my_imp_ind]; cerr << endl;}

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
#define rall(vec) vec.rbegin(), vec.rend()
#define all(vec) vec.begin(), vec.end()
#define eb emplace_back
#define pb push_back
#define pob pop_back
#define pf push_front
#define pof pop_front

std::mt19937_64 rng_64( std::chrono::steady_clock::now().time_since_epoch().count() );
constexpr int ilog2( int num ) { return num<=0?-1:8*sizeof(int) - __builtin_clz( num ) - 1; }
template<typename tpow,typename texp=int64_t> constexpr tpow mpow(tpow x,uint64_t e,tpow m){tpow res=1;while(e){if(e&1)res=(texp(res)*x)%m;e>>=1;x=(texp(x)*x)%m;}return res;}

/* Author: Oscar Vargas Pabon
Lo probe en https://codeforces.com/problemset/problem/2026/E 
Mi implementacion interpreta los valores de 'graph' como indices de 'edges'
Caso general funciona en $O(V*E*lgU)$, donde U es la capacidad maxima
    caso capacidad unitaria funciona en $O(min(E\sqrt{E},EV^{2/3}))$
    caso red unitaria (cada nodo se responsabiliza por 1 sola arista de
        capacidad unitaria) funciona en $O(E\sqrt{V})$
Es facil reducirlo a la version usual de dinitz

Remember,remember. Whenever we have a planar graph, a path of its dual graph represents
    bijectively a cut on the original graph. Essentially solves with dijkstra min-cut
Anfelesan taught me this one in regards to https://codeforces.com/gym/106178/problem/G
*/ template <typename ftype, ftype INF> struct Dinitz{
    constexpr static ftype inf(){return INF;}
    struct Edge { int u, v;
        ftype c1, c2, f;// c1 es la capacidad u->v y c2 es para v->u
        Edge()=default;
        Edge( int uu, int vv, ftype cc1, ftype cc2=0 ) : u(uu),v(vv),c1(cc1),c2(cc2),f(0) {};
        inline void push( int nd, ftype pushed ) { // empujar un flujo 'pushed' desde 'nd'
            if ( nd == u ) c1 -= pushed, c2 += pushed, f += pushed;
            else c1 += pushed, c2 -= pushed, f -= pushed;
        } inline int to( int nd ) const { // el opuesto a 'nd'
            return ( u == nd ) ? v : u;
        } inline ftype cap( int nd ) const { // la capacidad desde 'nd'
            return ( this->u == nd ) ? c1 : c2;
        } inline ftype flow( int nd ) const { // el flujo desde 'nd'
            return ( u == nd ) ? f : -f; }
    }; // level -> para el 'layered network'
    // blocked-> para omitir vertices 'blocked'
    // ptr -> para omitir aristas 'blocked'
    vector<int> level,ptr,q; vector<bool>blocked;

    vector<vector<int>> graph; vector<Edge> edges;
    int source, sink; Dinitz(int n){
        graph.resize(n); level.resize(n); q.resize(graph.size());
        ptr.resize(n); blocked.resize(n);
    } void add_edge( int u, int v, ftype c1, ftype c2=0 ) {
        // u->v has capacity c1; v->u has capacity c2   
        graph[u].push_back( edges.size() );
        graph[v].push_back( edges.size() );
        edges.emplace_back( Edge( u, v, c1, c2 ) );
    } ftype push_dfs( int nd, ftype flow=INF ) {
// Empuja el flujo por todas las aristas del 'layered graph' que pueda 
        if ( nd == sink || flow==0 ) return flow; // ya llegamos o no podemos empujar mas
        
        ftype push_flow = 0, edge_flow; int edge, nxt;
        while ( ptr[nd] < int(graph[nd].size()) && flow > push_flow ) {
            edge_flow = 0; edge = graph[nd][ptr[nd]]; nxt = edges[edge].to(nd);
            // si la arista pertenece al 'layered graph' y el vertice NO esta bloqueado
            if ( level[nd] < level[nxt] && !blocked[nxt] )  {
                edge_flow = push_dfs( nxt, min(edges[edge].cap(nd),flow-push_flow) );
                push_flow += edge_flow;
                edges[edge].push( nd, edge_flow );
            } ++ptr[nd];
        } --ptr[nd]; // para corregir el ultimo ++ptr[nd]
        blocked[nd] = (push_flow == 0); // verifica si se bloqueo el vertice
        return push_flow;
    } ftype max_flow( ) {
// Hace el maximo flujo del grafo (retorna el maxFlow, las asignaciones quedan en las aristas)
//      Funciona en peor caso $O(V*E*lgU)$
        ftype mx=0;for(const Edge&e:edges){//quitar esto para reducir
            if(mx<e.c1)mx=e.c1;if(mx<e.c2)mx=e.c2;// a la version
        }int ex=0;while((mx>>ex)>1)++ex;// usual de dinitz
        ftype flow=0;for(;ex>=0;--ex)do{
            fill(all(level),0);level[source]=1;
            int iq,nq=1;q[iq=0]=source;while(iq<nq&&!level[sink]){
                const int nd=q[iq++]; for(int eind:graph[nd]){
                    const int e=edges[eind].to(nd);
                    if(!level[e]&&(edges[eind].cap(nd)>>ex))
                        level[e]=level[nd]+1,q[nq++]=e;
                }
            } if(!level[sink])continue;
            fill(all(blocked),0); fill(all(ptr),0);
            flow+=push_dfs( source );
        } while(level[sink]); return flow; }
};

const int max_n = 1e6;
int b[max_n];
int s[max_n],t[max_n],l[max_n],u[max_n],c[max_n];

void solve() {
	int n,m;cin>>n>>m;
    rep(i,0,n)cin>>b[i];
    rep(i,0,m)cin>>s[i]>>t[i]>>l[i]>>u[i]>>c[i];
    Dinitz<int64_t,int64_t(1e18+33)> flow(n+2);
    flow.source=n;flow.sink=n+1;
    rep(i,0,n){
        if(b[i]>0)flow.add_edge(i,flow.sink,b[i]);
        else   flow.add_edge(flow.source,i,-b[i]);
    }

} int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;
#ifdef OSVARP
    auto start = chrono::high_resolution_clock::now();
#endif

    int t = 1;
    // cin >> t; 
    ++t; while ( --t ) solve();

#ifdef OSVARP
    auto end = chrono::high_resolution_clock::now();
    struct rusage usage; getrusage(RUSAGE_SELF, &usage);
    cerr << "\n<Execution time: "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms>\n<Memory used: "
        << usage.ru_maxrss << " kilobytes>" << endl;
#endif
    return 0; }