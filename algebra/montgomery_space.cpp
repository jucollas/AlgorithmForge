/*
Author: Oscar Vargas Pabon

It is though to work modulo primes, so inv (.inv,/,/=) and pow may not
	work properly otherwise.

I assume from my template :
inv :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}	

Tested in testing/test_alghelp.cpp and in fft/ntt stuff

Multiplications cost less (as we dont take costly modulus), but
	transforming to montgomery space now implies taking a modulus and
	detransforming from montgomery space implies the same work as a 
	multiplication
	
Tener cuidado si tengo que comparar con muchos enteros o cosas por el estilo

Adapted from https://en.algorithmica.org/hpc/number-theory/montgomery/
*/

typedef unsigned int uint; typedef unsigned long long ulint;
constexpr uint constexpr_calc_mr(uint n,int m_pow){uint nr=1;for(int i=0;i<ilog2(m_pow);++i)nr*=2-n*+nr;return nr;}
template < ulint raw_m, int m_pow=32, typename tint=uint,bool arbi_ntt=0 >
struct montgomery_int{
	constexpr static tint m=raw_m; static_assert(m>0);
	static_assert(m_pow>0&&m_pow<=32); 
	static constexpr tint mr=constexpr_calc_mr(m,m_pow);
	constexpr static tint mod(){return m;}
	constexpr static bool arbitrary_ntt(){return arbi_ntt;}
	
	constexpr static tint reduce(ulint x)noexcept{
		
		tint q= tint(x)*mr;
		tint y=((ulint)q*m)>>m_pow;
		
		tint res=(x>>m_pow)+m-y;
		if (res >= m) res -= m;
		return res;
	}
	constexpr static tint transform (tint x)noexcept{return (ulint(x)<<m_pow)%m;}
	constexpr static tint itransform(tint x)noexcept{return reduce(x);}
	
	tint vl;
	
	constexpr montgomery_int()noexcept:vl(0){};
	constexpr montgomery_int( int v)noexcept:vl(transform(v>=0?(v< int(m)?v:v%m):(m+v>=0?m+v:(v%m)+m))){};
	constexpr montgomery_int(lint v)noexcept:vl(transform(v>=0?(v<lint(m)?v:v%m):(m+v>=0?m+v:(v%m)+m))){};
	constexpr montgomery_int( uint v)noexcept:vl(transform(v<m?v:v%m)){};
	constexpr montgomery_int(ulint v)noexcept:vl(transform(v<ulint(m)?v:v%m)){};
	
	constexpr montgomery_int &operator +=(const montgomery_int &ot){ vl+=ot.vl; if(vl>=m)vl-=m; return *this; }
	constexpr montgomery_int  operator + (const montgomery_int &ot)const{ return montgomery_int(*this)+=ot; }
	constexpr montgomery_int &operator -=(const montgomery_int &ot){ vl=(vl>=ot.vl)?vl-ot.vl:vl+m-ot.vl; return *this; }
	constexpr montgomery_int  operator - (const montgomery_int &ot)const{ return montgomery_int(*this)-=ot; }
	constexpr montgomery_int &operator *=(const montgomery_int &ot){ vl= reduce((ulint)vl*ot.vl); return *this; }
	constexpr montgomery_int  operator * (const montgomery_int &ot)const{ return montgomery_int(*this)*=ot; }
	constexpr montgomery_int &operator /=(const montgomery_int &ot){ (*this)*=ot.inv(); return *this; }
	constexpr montgomery_int  operator / (const montgomery_int &ot)const{ return montgomery_int(*this)/=ot; }
	
	constexpr montgomery_int operator -()const {return montgomery_int(0)-(*this);}
	constexpr montgomery_int pow(lint e)const{
		montgomery_int rs=1,ac=*this;while(e){
			if(e&1ll)rs*=ac;
			e>>=1;ac*=ac;
		}
		return rs;
	}
	constexpr montgomery_int inv()const{return this->pow(m-2);}//Fermats little theorem
	
	constexpr bool operator ==(const montgomery_int &ot)const{return vl==ot.vl;} 
	constexpr bool operator !=(const montgomery_int &ot)const{return vl!=ot.vl;}
	constexpr bool operator ==(const  uint &ot)const{return reduce(vl)==ot;}
	constexpr bool operator ==(const  int &ot )const{return reduce(vl)==ot;}
	
	constexpr operator bool() const { return itransform(vl); }
	constexpr operator  int() const { return itransform(vl); }
	constexpr operator long long() const{ return int(*this); }
	
	friend ostream &operator<<(ostream &os,const montgomery_int &ac){return os << itransform(ac.vl);}
	friend istream &operator>>(istream&is,montgomery_int &ac){int v;is>>v;ac=montgomery_int(v);return is;}	
};  typedef montgomery_int<mod> mint;

// for a more focalized use
struct montgomery_space{
	static const int m_pow=32;static_assert(m_pow>0&&m_pow<=32);
	uint m,mr;
	constexpr montgomery_space(uint md):m(md),mr(1){
		mr=constexpr_calc_mr(m,m_pow);
	}
	
	uint reduce(ulint x)const{
		
		uint q= uint(x)*mr;
		uint y=((ulint)q*m)>>m_pow;
		
		uint res=(x>>m_pow)+m-y;
		if (res >= m) res -= m;
		return res;
	}
	uint transform (uint x)const{return (ulint(x)<<m_pow)%m;}
	uint itransform(uint x)const{return reduce(x);}
	
	uint mul(uint x,uint y)const {return reduce((ulint)x*y);}
	uint add(uint x,uint y)const {uint rs=x+y;if(rs>=m)rs-=m;return rs;}
	uint sub(uint x,uint y)const {return (x>=y)?x-y:x+m-y;}
	
};