/*
Author: Oscar Vargas Pabon

It is though to work modulo primes, so inv (.inv,/,/=) and pow may not
	work properly otherwise.

I assume from my template :
inv :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}	

Tested in testing/test_alghelp.cpp and in fft/ntt stuff
*/

template<lint raw_m,typename tint=int,bool arbi_ntt=0>
struct modulo_int{
	constexpr static tint m=raw_m; static_assert(m>0);
	constexpr static tint mod(){return m;}
	constexpr static bool arbitrary_ntt(){return arbi_ntt;}
	
	tint vl;
	constexpr modulo_int()noexcept:vl(0){};
	constexpr modulo_int( int v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(v%m)+m)){};
	constexpr modulo_int(lint v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(v%m)+m)){};
	constexpr modulo_int(unsigned long long v)noexcept:vl(v<m?v:v%m){};
	
	constexpr modulo_int &operator +=(const modulo_int &ot){ vl+=ot.vl; if(vl>=m)vl-=m; return *this; }
	constexpr modulo_int  operator + (const modulo_int &ot)const{ return modulo_int(*this)+=ot; }
	constexpr modulo_int &operator -=(const modulo_int &ot){ vl=(vl>ot.vl)?vl-ot.vl:vl+m-ot.vl; return *this; }
	constexpr modulo_int  operator - (const modulo_int &ot)const{ return modulo_int(*this)-=ot; }
	constexpr modulo_int &operator *=(const modulo_int &ot){ vl=(vl*1ll*ot.vl)%m; return *this; }
	constexpr modulo_int  operator * (const modulo_int &ot)const{ return modulo_int(*this)*=ot; }
	constexpr modulo_int &operator /=(const modulo_int &ot){ (*this)*=ot.inv(); return *this; }
	constexpr modulo_int  operator / (const modulo_int &ot)const{ return modulo_int(*this)/=ot; }
	
	constexpr modulo_int inv()const{return modulo_int(vl).pow(m-2);}//Fermats little theorem
	constexpr modulo_int operator -()const {return modulo_int(-vl);}
	constexpr modulo_int pow(lint e)const{return modulo_int(mpow<tint>(vl,e%(m-1),m));}
	
	constexpr bool operator ==(const modulo_int &ot)const{return vl==ot.vl;} 
	constexpr bool operator ==(const  int &ot)const{return vl==ot;   }
	constexpr bool operator !=(const modulo_int &ot)const{return vl!=ot.vl;}
	
	constexpr operator bool() const { return vl; }
	constexpr operator  int() const { return vl; }
	constexpr operator long long() const { return vl; }
	
	friend ostream &operator<<(ostream &os,const modulo_int &ac){return os << ac.vl;}
	friend istream &operator>>(istream&is,modulo_int &ac){int v;is>>v;ac=modulo_int(v);return is;}	
}; typedef modulo_int<mod> mint;
