/*
Author: Oscar Vargas Pabon

It is though to work modulo primes, so inv (.inv,/,/=) and pow may not
	work properly otherwise.

I assume from my template :
inv :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}	

Tested in testing/test_alghelp.cpp and in fft/ntt stuff
*/
template<unsigned long long raw_m,typename tint=unsigned int,typename tmul=unsigned long long,bool arbi_ntt=0>
struct modulo_int{ constexpr static tint m=raw_m; static_assert(m>0);
	constexpr static tint mod(){return m;}
	constexpr static bool arbitrary_ntt(){return arbi_ntt;}
	
	tint vl;
	inline constexpr modulo_int()noexcept:vl(0){};
	inline constexpr modulo_int(      int v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(v%m)+m)){};
	inline constexpr modulo_int(long long v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(v%m)+m)){};
	inline constexpr modulo_int(unsigned       int v)noexcept:vl(v<m?v:v%m){};
	inline constexpr modulo_int(unsigned long long v)noexcept:vl(v<m?v:v%m){};
	
	inline constexpr modulo_int &operator +=(const modulo_int &ot){ vl= m-vl>ot.vl?vl+ot.vl:ot.vl-(m-vl); return *this; }
	inline constexpr modulo_int  operator + (const modulo_int &ot)const{ return modulo_int(*this)+=ot; }
	inline constexpr modulo_int &operator -=(const modulo_int &ot){ vl=(vl>=ot.vl)?vl-ot.vl:(m-ot.vl)+vl; return *this; }
	inline constexpr modulo_int  operator - (const modulo_int &ot)const{ return modulo_int(*this)-=ot; }
	inline constexpr modulo_int &operator *=(const modulo_int &ot){ vl=(tmul(vl)*ot.vl)%m; return *this; }
	inline constexpr modulo_int  operator * (const modulo_int &ot)const{ return modulo_int(*this)*=ot; }
	inline constexpr modulo_int &operator /=(const modulo_int &ot){ (*this)*=ot.inv(); return *this; }
	inline constexpr modulo_int  operator / (const modulo_int &ot)const{ return modulo_int(*this)*=ot.inv(); }
	
	inline constexpr modulo_int operator -()const {return modulo_int(-vl);}
	inline constexpr modulo_int inv()const{return mpow<tint,tmul>(vl,m-2,m);}//Fermats little theorem
	inline constexpr modulo_int pow(long long e)const{return e>=0?modulo_int(mpow<tint,tmul>(vl,e%(m-1),m)):inv().pow(-e);}
	
	inline constexpr bool operator ==(const modulo_int &ot)const{return vl==ot.vl;} 
	inline constexpr bool operator ==(const  int &ot)const{return vl==ot;   }
	inline constexpr bool operator !=(const modulo_int &ot)const{return vl!=ot.vl;}
	
	inline constexpr operator bool() const { return vl; }
	inline constexpr operator  int() const { return vl; }
	inline constexpr operator unsigned int() const{return vl;}
	inline constexpr operator long long() const { return vl; }
	inline constexpr operator unsigned long long() const{return vl;}
	
	friend ostream &operator<<(ostream &os,const modulo_int &ac){return os << ac.vl;}
	friend istream &operator>>(istream&is,modulo_int &ac){int v;is>>v;ac=modulo_int(v);return is;}	
}; typedef modulo_int<mod> mint;
