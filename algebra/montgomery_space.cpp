/*
Autor: Oscar Vargas Pabon

Las multiplicaciones son menos costosas, pero transformar y detransformar del
	espacio montgomery cuesta relativamente bastante.
	
	Tener cuidado si tengo que comparar con muchos enteros o cosas por el estilo

Adaptado de https://en.algorithmica.org/hpc/number-theory/montgomery/
*/

typedef unsigned int uint; typedef unsigned long long ulint;
constexpr uint constexpr_calc_modr(uint n,int m_pow){uint nr=1;for(int i=0;i<ilog2(m_pow);++i)nr*=2-n*+nr;return nr;}
struct montgomery_int{
	static const int m_pow=32;static_assert(m_pow>0&&m_pow<=32);
	static const uint mod= 998244353; static_assert(mod>0);
	static constexpr uint modr=constexpr_calc_modr(mod,m_pow);
	
	constexpr static uint reduce(ulint x)noexcept{
		
		uint q= uint(x)*modr;
		uint m=((ulint)q*mod)>>m_pow;
		
		uint res=(x>>m_pow)+mod-m;
		if (res >= mod) res -= mod;
		return res;
	}
	constexpr static uint transform (uint x)noexcept{return (ulint(x)<<m_pow)%mod;}
	constexpr static uint itransform(uint x)noexcept{return reduce(x);}
	
	uint vl;
	
	constexpr montgomery_int()noexcept:vl(0){};
	constexpr montgomery_int( int v)noexcept:vl(transform(v>=0?v%mod:(v%mod)+mod)){};
	constexpr montgomery_int(lint v)noexcept:vl(transform(v>=0?v%mod:(v%mod)+mod)){};
	constexpr montgomery_int(uint v)noexcept:vl(transform(v%mod)){};
	constexpr montgomery_int(ulint v)noexcept:vl(transform(v%mod)){};
	
	montgomery_int &operator +=(const montgomery_int &ot){ vl+=ot.vl; if(vl>=mod)vl-=mod; return *this; }
	montgomery_int  operator + (const montgomery_int &ot)const{ return montgomery_int(*this)+=ot; }
	montgomery_int &operator -=(const montgomery_int &ot){ vl=(vl>=ot.vl)?vl-ot.vl:vl+mod-ot.vl; return *this; }
	montgomery_int  operator - (const montgomery_int &ot)const{ return montgomery_int(*this)-=ot; }
	montgomery_int &operator *=(const montgomery_int &ot){ vl= reduce((ulint)vl*ot.vl); return *this; }
	montgomery_int  operator * (const montgomery_int &ot)const{ return montgomery_int(*this)*=ot; }
	montgomery_int &operator /=(const montgomery_int &ot){ (*this)*=ot.inverse(); return *this; }
	montgomery_int  operator / (const montgomery_int &ot)const{ return montgomery_int(*this)/=ot; }
	
	montgomery_int operator -()const {return montgomery_int(0)-(*this);}
	montgomery_int pow(lint e)const{
		montgomery_int rs=1,ac=*this;while(e){
			if(e&1ll)rs*=ac;
			e>>=1;ac*=ac;
		}
		return rs;
	}
	montgomery_int inverse()const{return this->pow(mod-2);}//Fermats little theorem
	
	bool operator ==(const montgomery_int &ot)const{return vl==ot.vl;} 
	bool operator ==(const  uint &ot)const{return vl==transform(ot%mod);}
	bool operator ==(const  int &ot )const{return vl==transform(ot%mod);}
	bool operator !=(const montgomery_int &ot)const{return vl!=ot.vl;}
	
	operator bool() const { return itransform(vl); }
	operator  int() const { return itransform(vl); }
	
	friend ostream &operator<<(ostream &os,const montgomery_int &ac){return os << itransform(ac.vl);}
	friend istream &operator>>(istream&is,montgomery_int &ac){int v;cin>>v;ac=montgomery_int(v);return is;}	
}; const int mod=montgomery_int::mod; typedef montgomery_int mint;

// for a more focalized use
struct montgomery_space{
	static const int m_pow=32;static_assert(m_pow>0&&m_pow<=32);
	uint mod,modr;
	constexpr montgomery_space(uint md):mod(md),modr(1){
		modr=constexpr_calc_modr(mod,m_pow);
	}
	
	uint reduce(ulint x)const{
		
		uint q= uint(x)*modr;
		uint m=((ulint)q*mod)>>m_pow;
		
		uint res=(x>>m_pow)+mod-m;
		// debug(x,"enestaver",q,m,res);
		if (res >= mod) res -= mod;
		return res;
	}
	uint transform (uint x)const{return (ulint(x)<<m_pow)%mod;}
	uint itransform(uint x)const{return reduce(x);}
	
	uint mul(uint x,uint y)const {return reduce((ulint)x*y);}
	uint add(uint x,uint y)const {uint rs=x+y;if(rs>=mod)rs-=mod;return rs;}
	uint sub(uint x,uint y)const {return (x>=y)?x-y:x+mod-y;}
	
};