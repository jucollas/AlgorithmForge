/*
Autor:Oscar Vargas Pabon

Modulo Integer
makes stuff a bit easier
*/

struct mint{
	const static int mod= 998244353; static_assert(mod>0);
	
	int vl;
	constexpr mint()noexcept:vl(0){};
	constexpr mint( int v)noexcept:vl(v>=0?v%mod:(v%mod)+mod){};
	constexpr mint(lint v)noexcept:vl(v>=0?v%mod:(v%mod)+mod){};
	constexpr mint(unsigned long long v)noexcept:vl(v%mod){};
	
	mint &operator +=(const mint &ot){ vl+=ot.vl; if(vl>=mod)vl-=mod; return *this; }
	mint  operator + (const mint &ot)const{ return mint(*this)+=ot; }
	mint &operator -=(const mint &ot){ vl-=ot.vl; if(vl<0)vl+=mod; return *this; }
	mint  operator - (const mint &ot)const{ return mint(*this)-=ot; }
	mint &operator *=(const mint &ot){ vl=(vl*1ll*ot.vl)%mod; return *this; }
	mint  operator * (const mint &ot)const{ return mint(*this)*=ot; }
	mint &operator /=(const mint &ot){ (*this)*=ot.inverse(); return *this; }
	mint  operator / (const mint &ot)const{ return mint(*this)/=ot; }
	
	mint inverse()const{return mint(mpow(vl,mod-2,mod));}//Fermats little theorem
	mint operator -()const {return mint(-vl);}
	mint pow(lint e)const{return mint(mpow(vl,e%(mod-1),mod));}
	
	bool operator ==(const mint &ot)const{return vl==ot.vl;} 
	bool operator ==(const  int &ot)const{return vl==ot;   }
	bool operator !=(const mint &ot)const{return vl!=ot.vl;}
	
	operator bool() const { return vl; }
	operator  int() const { return vl; }
	
	friend ostream &operator<<(ostream &os,const mint &ac){return os << ac.vl;}
	friend istream &operator>>(istream&is,mint &ac){int v;cin>>v;ac=mint(v);return is;}	
}; const int mod=mint::mod;