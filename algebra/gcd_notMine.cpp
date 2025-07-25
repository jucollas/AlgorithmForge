// https://cp-algorithms.com/algebra/extended-euclid-algorithm.html
int gcd(int a, int b, int& x, int& y) {
	// a*x + b*y = gcd(a,b)
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int x1, y1;
    int d = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}