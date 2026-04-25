/*
Autor: Oscar vargas Pabon

Nota: en vez de precomputar ebase podriSa usar mpow para hacerlo en O(log n)

NO OLVIDAR INVOCAR 'precompute_prime' ANTES DE USAR

Testeado en; https://codeforces.com/problemset/problem/2132/G
			https://codeforces.com/problemset/problem/727/E
*/

// const int template_limit=1e6;
// mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
// #define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
// template<typename tpow> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}
// asumo el anterior de mi template

const int N_HASH=2, PRIME_TESTING=3;
const int LH_B=1e9,RH_B=1<<31;
int base[N_HASH],modul[N_HASH],ebase[template_limit][N_HASH];

void precompute_prime(){
  auto rgen = [&](int lb, int rb ) -> int {
    return lb+(rng_64())%(rb-lb+1);
  }; for ( int i = 0 ; i < N_HASH ; ++i ) base[i]=rgen(2,LH_B);
  int posi; bool prim;
  for ( int i = 0 ; i < N_HASH ; ++i ) {
    do {
      posi = rgen(LH_B,RH_B);
      prim = mpow<int>(2,posi-1,posi)==1;
      for (int j = 0 ; j < N_HASH && prim ; ++j )
        prim=mpow<int>(base[i],posi-1,posi)==1;
      for (int j = 0 ; j < PRIME_TESTING && prim ; ++j )
        prim=mpow<int>(rgen(2,LH_B),posi-1,posi)==1;
    } while ( !prim ) ;
    modul[i] = posi;
  } // precomputar las potencias de la base _ tambien podria hacerlo en O(log n) con mpow
  rep(j,0,N_HASH) ebase[0][j]=1;
  rep(i,0,template_limit-1) rep(j,0,N_HASH) ebase[i+1][j]=(ebase[i][j]*1ll*base[j])%modul[j];
}

struct Hash{
  int len = 0; array<int,N_HASH> h;
  Hash() :len(0) {rep(i,0,N_HASH)h[i]=0;};
  Hash(const string &cad) { len = cad.size();
    for ( int i = 0 ; i < N_HASH ; ++i ) { h[i]=0;
      for ( char c : cad ) h[i] = ( h[i]*1ll*base[i] + int(c) )%modul[i];
    }
  } Hash( char c ) { len=1;
	rep(i,0,N_HASH)h[i]= int(c)%modul[i];  
  } bool operator == ( const Hash &o ) const { return len==o.len&&h==o.h; }
  bool operator < ( const Hash &o ) const { return h<o.h; }
  Hash operator + ( const Hash &o ) const { // *this + o
    Hash res; for ( int i = 0 ; i < N_HASH ; ++i ) {
      res.h[i] = ( ebase[o.len][i]*1ll*h[i] + o.h[i] )%modul[i];
    } res.len = len+o.len; return res;
  } void lconc( char c ) { // hago this = c + this
	rep(i,0,N_HASH) h[i] = ( ebase[len][i]*1ll*lint(c) + h[i] ) % modul[i];
	++len;
  } void rconc( char c ) { // hago this=this+c
	rep(i,0,N_HASH) h[i] = ( h[i]*1ll*base[i] + lint(c))%modul[i];
	++len;
  } void ldeconc ( const Hash &l ) { // quito el prefijo l de *this
    for ( int i = 0 ; i < N_HASH ; ++i ) {
      h[i] -=(l.h[i]*1ll*ebase[len-l.len][i] )%modul[i];
	  if(h[i]<0)h[i]+=modul[i];
    } len -= l.len;
  } void rdeconc( const Hash &r ) { // quito el sufijo r de *this
    for ( int i = 0 ; i < N_HASH ; ++i ) {
      h[i] = ( (h[i]-r.h[i])*1ll*mpow<int>(ebase[r.len][i], modul[i]-2,modul[i]) ) %modul[i];
	  if(h[i]<0)h[i]+=modul[i];
    } len -=r.len;
  }
};