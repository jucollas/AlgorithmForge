/*
Autor:Oscar Vargas Pabon

Modulo Integer
makes stuff a bit easier
*/
struct modulo_int{
	static const int mod= 998244353; static_assert(mod>0);
	
	int vl;
	constexpr modulo_int()noexcept:vl(0){};
	constexpr modulo_int( int v)noexcept:vl(v>=0?v%mod:(v%mod)+mod){};
	constexpr modulo_int(lint v)noexcept:vl(v>=0?v%mod:(v%mod)+mod){};
	constexpr modulo_int(unsigned long long v)noexcept:vl(v%mod){};
	
	modulo_int &operator +=(const modulo_int &ot){ vl+=ot.vl; if(vl>=mod)vl-=mod; return *this; }
	modulo_int  operator + (const modulo_int &ot)const{ return modulo_int(*this)+=ot; }
	modulo_int &operator -=(const modulo_int &ot){ vl-=ot.vl; if(vl<0)vl+=mod; return *this; }
	modulo_int  operator - (const modulo_int &ot)const{ return modulo_int(*this)-=ot; }
	modulo_int &operator *=(const modulo_int &ot){ vl=(vl*1ll*ot.vl)%mod; return *this; }
	modulo_int  operator * (const modulo_int &ot)const{ return modulo_int(*this)*=ot; }
	modulo_int &operator /=(const modulo_int &ot){ (*this)*=ot.inverse(); return *this; }
	modulo_int  operator / (const modulo_int &ot)const{ return modulo_int(*this)/=ot; }
	
	modulo_int inverse()const{return modulo_int(mpow(vl,mod-2,mod));}//Fermats little theorem
	modulo_int operator -()const {return modulo_int(-vl);}
	modulo_int pow(lint e)const{return modulo_int(mpow(vl,e%(mod-1),mod));}
	
	bool operator ==(const modulo_int &ot)const{return vl==ot.vl;} 
	bool operator ==(const  int &ot)const{return vl==ot;   }
	bool operator !=(const modulo_int &ot)const{return vl!=ot.vl;}
	
	operator bool() const { return vl; }
	operator  int() const { return vl; }
	
	friend ostream &operator<<(ostream &os,const modulo_int &ac){return os << ac.vl;}
	friend istream &operator>>(istream&is,modulo_int &ac){int v;cin>>v;ac=modulo_int(v);return is;}	
}; const int mod=modulo_int::mod; typedef modulo_int mint;