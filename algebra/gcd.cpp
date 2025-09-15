/*
Plagiado epicamente de cp-algorithms

se me olvido como funciona la extendida de euclides.
------------------UNFINISHED
g:=gcd(a,b); WLOG a>=b; a=(a/b)+(a%b); then

	g=x*a+y*b
	g=x*((a/b)+(a%b))+y*b
	g=
------------------

*/

lint gcd(lint a, lint b, lint& x, lint& y) {
	// a*x+b*y==gcd(a,b);
    if (b == 0ll) {
        x = 1;
        y = 0;
        return a;
    }
    lint x1, y1;
    lint d = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}