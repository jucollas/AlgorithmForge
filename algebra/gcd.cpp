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
template<typename tgcd>
tgcd gcd(tgcd a, tgcd b, tgcd& x, tgcd& y) {
	// a*x+b*y==gcd(a,b);
    if (b == tgcd(0)) {
        x = 1;
        y = 0;
        return a;
    }
    tgcd x1, y1;
    tgcd d = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}