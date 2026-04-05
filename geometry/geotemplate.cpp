/*
Autor: Oscar Vargas Pabon

GRAFOS PLANARES:::::
c:componentes conexos ; R:regiones
m:aristas; n:nodos
R+n=m+c+1

Probado en UVA 13117 y 11894
*/
const double eps=1e-8;
template<typename Tpt> struct Point{
	Tpt x,y;
	Point()=default;
	Point(Tpt x, Tpt y ) : x(x),y(y){};
	Point rot90(){return{-y,x};}
	
	// NOTA: puedo querer prescindir del sqrt
	Tpt norm()const{return sqrt(x*x+y*y);}
	Tpt norm2()const{return x*x+y*y;}
	
	Point scale(Tpt v)const{return{v*x,v*y};}
	Point operator -()const{return{-x,-y};}
	Point operator -(const Point&o)const{return{x-o.x,y-o.y};}
	// suma de vectores
	Point operator +(const Point&o)const{return{x+o.x,y+o.y};}
	
	// point product; 0->ortogonal; +->same dir ; - ->opposite dir
	Tpt operator *(const Point&o)const{return x*o.x+y*o.y;}
	// cross product; 0->colinear; +->left side; - ->right side
	Tpt operator %(const Point&o)const{return x*o.y-y*o.x;}

// abs(v%u)==u.norm2()*v.norm2()*sin(theta)
// (v*u)==u.norm2()*v.norm2()*cos(theta)
// v%u == area paralelogramo delimitado por u y v
	bool operator<(const Point &o)const{return make_pair(x,y)<make_pair(o.x,o.y);}
	bool operator==(const Point&o)const{return x==o.x&&y==o.y;}
	bool operator!=(const Point&o)const{return!(*this==o);}
}; typedef Point<double> Pt;
ostream & operator << (ostream &out, const Pt &p){ out << "("<<p.x<<","<<p.y<<")"; return out; }

template<typename Tpt> struct Line{
	Point<Tpt> p,v;// L(t)=p+tv
	Line()=default;
	
	// L(0)=a; L(1)=b; L(t) 0<=t<=1 esta en el segmento ab
	Line(Point<Tpt> a,Point<Tpt> b):p(a),v(b-a){};
	Point<Tpt> eval(Tpt t) { return p+v.scale(t);}
	// hallo t tal que L(t)== interseccion entre this y o
	// recordar que puedo trabajar en enteros hasta la division
	Tpt inter(const Line&o){return{((p-o.p)%o.v)/(v%o.v)};}
	
	// dados v,u vectores; el valor h donde {0,v(h),u} forman un triangulo rectangulo
	// con con angulo recto A{o,v(h),u} =90 es tal que u*v=h*(v*v)
};typename Line<double> Ln;

ostream & operator << (ostream &out, const Ln &l){ out << "<L(t)="<<l.p<<"+t"<<l.v<<">"; return out; }