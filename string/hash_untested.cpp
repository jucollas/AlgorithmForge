/*
Autor: Oscar vargas Pabon

Nota: en vez de precomputar ebase podría usar mpow para hacerlo en O(log n)
*/


lint mpow( lint x, lint e, lint m ) {
  lint res = 1ll;
  while ( e ) {
    if ( e&1ll ) res = (res*x)%m;
    e >>= 1; x = (x*x)%m;
  }
  return res;
}

//defined in header <random> :- mersenne-twister
mt19937_64 random_64( chrono::steady_clock::now().time_since_epoch().count() );
// #define rep(i,n) for(int i = 0 ; i < n ; ++i )
// asumo el anterior de mi template

const int N_HASH=2, PRIME_TESTING=3;
const lint LH_B=1e9,RH_B=1ll<<31;
lint base[N_HASH],modul[N_HASH],ebase[template_limit][N_HASH];

void precompute_prime(){
  function<lint(lint,lint)> rgen = [&](lint lb, lint rb ) {
    return lb+(random_64())%(rb-lb+1ll);
  };
  for ( int i = 0 ; i < N_HASH ; ++i ) base[i]=rgen(2,LH_B);
  lint posi; bool prim;
  for ( int i = 0 ; i < N_HASH ; ++i ) {
    do {
      posi = rgen(LH_B,RH_B);
      prim = mpow(2,posi-1,posi)==1ll;
      for (int j = 0 ; j < N_HASH && prim ; ++j )
        prim=mpow(base[i],posi-1,posi)==1ll;
      for (int j = 0 ; j < PRIME_TESTING && prim ; ++j )
        prim=mpow(rgen(2,LH_B),posi-1,posi)==1ll;
    } while ( !prim ) ;
    modul[i] = posi;
  }
  // precomputar las potencias de la base _ tambien podria hacerlo en O(log n) con mpow
  rep(j,N_HASH) ebase[0][j]=1;
  rep(i,template_limit-1) rep(j,N_HASH) ebase[i+1][j]=(ebase[i][j]*base[j])%modul[j];
}

class Hash{
public:
  array<lint,N_HASH> h;
  int len = 0;
  Hash() :len(0),h({0,0}) {};
  Hash(const string &cad) {
    len = cad.size();
    for ( int i = 0 ; i < N_HASH ; ++i ) {
      h[i]=0ll;
      for ( char c : cad ) h[i] = ((h[i]*base[i])%modul[i] + lint(c) )%modul[i];
    }
  }
  Hash( char c ) {
	len=1;
	rep(i,N_HASH)h[i]= lint(c)%modul[i];  
  }
  bool operator == ( const Hash &o ) const { return len==o.len&&h==o.h; }
  Hash operator + ( const Hash &o ) const {
    // *this + o
    Hash res;
    for ( int i = 0 ; i < N_HASH ; ++i ) {
      res.h[i] = ( (ebase[o.len][i]*h[i])%modul[i] + o.h[i] )%modul[i];
    }
    res.len = len+o.len;
    return res;
  }
  void lconc( char c ) {
	// hago this = c + this
	rep(i,N_HASH) h[i] = ( ebase[len][i]*lint(c) + h[i] ) % modul[i];
	++len;
  }
  void rconc( char c ) {
	  // hago this=this+c
	rep(i,N_HASH) h[i] = ( h[i]*base[i] + lint(c))%modul[i];
	++len;
  }
  void ldeconc ( const Hash &l ) {
    // quito el prefijo l de *this
	// cout << l.len << " _ " << endl;
    for ( int i = 0 ; i < N_HASH ; ++i ) {
      h[i] -=(l.h[i]*ebase[len-l.len][i] )%modul[i];
	  h[i] = ( h[i] + modul[i] ) %modul[i];
    }
    len -= l.len;
  }
  void rdeconc( const Hash &r ) {
    // quito el sufijo r de *this
    for ( int i = 0 ; i < N_HASH ; ++i ) {
      h[i] = ( (h[i]-r.h[i]) * mpow(ebase[r.len][i], modul[i]-2,modul[i]) ) %modul[i];
	  h[i] = ( h[i] + modul[i] ) %modul[i];
    }
    len -=r.len;
  }
};