// Ill eventually try to understand this shit
///// I think it is the key to aladdin from COCI 2009/2010 that I found in the CCPL
/// floor sum shenanigans



// taken from https://asfjwd.github.io/2020-04-24-floor-sum-ap/
// its not safe for cases when a<0 or b<0 xdxdxd
long long FloorSumAP(long long a, long long b, long long c, long long n){}
// calculating $\sum_{x=0}^n \lfloor (ax+b)/c \rfloor$
  if(!a) return (b / c) * (n + 1);
  if(a >= c || b >= c) return ( ( n * (n + 1) ) / 2) * (a / c) + (n + 1) * (b / c) + FloorSumAP(a % c, b % c, c, n);
  long long m = (a * n + b) / c;
  return m * n - FloorSumAP(c, c - b - 1, a, m - 1);
}
// taken from  atcoder library
namespace internal{
constexpr long long safe_mod(long long x,long long m){ x%=m;if(x<0)x+=m;return x; }
unsigned long long floor_sum_unsigned(unsigned long long n,
                                      unsigned long long m,
                                      unsigned long long a,
                                      unsigned long long b) {
    unsigned long long ans = 0;
    while (true) {
        if (a >= m) {
            ans += n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if (b >= m) {
            ans += n * (b / m);
            b %= m;
        }

        unsigned long long y_max = a * n + b;
        if (y_max < m) break;
        // y_max < m * (n + 1)
        // floor(y_max / m) <= n
        n = (unsigned long long)(y_max / m);
        b = (unsigned long long)(y_max % m);
        std::swap(m, a);
    }
    return ans;
}

}//namespace internal
long long floor_sum(long long n, long long m, long long a, long long b) {
	// Computes $\sum_{i=0}^{n-1} \lfloor\frac{a*i+b}{m}\rfloor$ in time O(lgn)
    assert(0 <= n && n < (1LL << 32));
    assert(1 <= m && m < (1LL << 32));
    unsigned long long ans = 0;
    if (a < 0) {
        unsigned long long a2 = internal::safe_mod(a, m);
        ans -= 1ULL * n * (n - 1) / 2 * ((a2 - a) / m);
        a = a2;
    }
    if (b < 0) {
        unsigned long long b2 = internal::safe_mod(b, m);
        ans -= 1ULL * n * ((b2 - b) / m);
        b = b2;
    }
    return ans + internal::floor_sum_unsigned(n, m, a, b);
}