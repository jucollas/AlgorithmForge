/*
Autor: Oscar Vargas Pabon

GRAFOS PLANARES:::::
c:componentes conexos ; R:regiones
m:aristas; n:nodos
R+n=m+c+1
*/
const double eps=1e-8;
typedef double Tpt;
class Pt{
	public:
	Tpt x,y;
	Pt()=default;
	Pt(Tpt x, Tpty ) : x(x),y(y){};
	Pt rot90(){return{-y,x};}
	// NOTA: puedo querer prescindir del sqrt
	Tpt norm2()const{return sqrt(x*x+y*y);}
	
	Pt scale(Tpt v)const{return{v*x,v*y};}
	Pt operator -()const{return{-x,-y};}
	Pt operator -(const Pt&o)const{return{x-o.x,y-o.y};}
	// suma de vectores
	Pt operator +(const Pt&o)const{return{x+o.x,y+o.y};}
	
	// point product; 0->ortogonal; +->same dir ; - ->opposite dir
	Tpt operator *(const Pt&o)const{return x*o.x+y*o.y;}
	// cross product; 0->colinear; +->left side; - ->right side
	Tpt operator %(const Pt&o)const{return x*o.y-y*o.x;}

// abs(v%u)==u.norm2()*v.norm2()*sin(theta)
// (v*u)==u.norm2()*v.norm2()*cos(theta)
// v%u == area paralelogramo delimitado por u y v
};

class Ln{
	public:
	Pt p,v;// L(t)=p+tv
	Ln()=default;
	
	// L(0)=a; L(1)=b; L(t) 0<=t<=1 esta en el segmento ab
	Ln(Pt a,Pt b):p(a),v(b-a){};
	Pt eval(Tpt t) { return p+v.scale(t);}
	// hallo t tal que L(t)== interseccion entre this y o
	// recordar que puedo trabajar en enteros hasta la division
	Tpt inter(const Ln&o){return{((p-o.p)%o.v)/(v%o.v)};}
};