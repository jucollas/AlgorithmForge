/*
Autor: Oscar Vargas Pabon
These are various tricks to be used at will
*/

// iterating over ranges [l_i,r_i]\subseteq[L,R] such that (assuming integer division)
//////   $ \forall_{j,h\in[l_i,r_i]}n/j=n/h $
/// can be seen in https://codeforces.com/problemset/problem/2072/G
void diviter(int n, int L, int R ){
	int l = L,r;
	while( l <= R ){
		int nh = n/l;
		r = min(R,n/nh);
		
		/// do stuff with the range $i\in[l,r] n/i=nh$
		
		l=r+1;
	}
}

// ham-path
// $phi(nd,msk)= strt\in N(nd)$ OR $OR_{e\in N(nd)} phi(nd,msk|(1<<nd))$
// In cases where I'm only looking for cycles containing msk, then I will start
// from every node; however, I don't need to erase the memory. The reason is that
// if such cycle exists I can find it with any $strt\in msk$, if there is none it won't
// matter which $strt\in msk$ I choose.
// This idea can be seen in https://codeforces.com/problemset/problem/1804/E


/////// XOR weird properties
// https://codeforces.com/blog/entry/141285
// a+b= a^b + 2(a&b)
// a^b <= 2*max(a,b)
// a^b=(a|b)-(a&b)
// a-b <= a^b <= a+b
// A[0..n) such that A[i]<=A[i+1], then min_{i!=j}(A[i]^A[j])=min_i(A[i]^A[i+1])

/////// Range Xor Shenanigans (assuming range is [l,r])
// Question: Is subset in range interesting? Is it doable in O(1)?
// I used this in industrial nim and in https://codeforces.com/problemset/problem/2056/F2

// calculating xor in range
int xor_range(int m){
	//calculates $XOR_{x=0}^m x$ in O(1). Observe it includes m.
	int res=(m&1)?0:m;
	return res^=((m+1)>>1)&1;
}

// calculating xor in range V.2
template<typename tXor> tXor xor_range_V2(tXor n){
	vector<tXor> _p = {n, 1, n + 1, 0}; return _p[n % 4];
}// This one follows by induction or showing $\Xor_{i=0}^{4n-1}i=0$
// Close inspection shows both versions do the same
// https://codeforces.com/blog/entry/141285
 
// The idea is that if I fix the bits that are always set, then I can exclude them
// x= 0010100
// y= ??1?1?? 
// meaning that if I exclude the 1's it becomes a xor_range problem
// this involves calculating the new m (since the actual m may not me superset)
//    and re-adding again the 'erased' bits after xor_range
int super_set(int m,int x){
	//calculates $XOR_{y=0}^m (y\&x==x)*y$ in O(lgU)
	if(x>m)return 0;
	
	//Finding last superset in range
	int lgi=ilog2(m-1),lst=x;
	rep(e,lgi,-1)if((lst|(1<<e))<m)lst|=1<<e;
	
	
	int neo=0,xx=x,e=0;
	while(lst){//Reducing to xor_range
		if(!(xx&1))neo|=(lst&1)<<e++;
		lst>>=1;xx>>=1;
	}
	int skw=xor_range(neo);xx=x;
	
	int res=0; e=0;
	while(skw){//reconstruct answer from xor_range
		if(!(xx&1)){
			res|=(skw&1)<<e;
			skw>>=1;
		}
		xx>>=1; ++e;
	}
	
	// add the bits in x (which are always set)
	if(x)res|=((neo+1)&1)*x;
	return res;
}

///// I think it is the key to aladdin from COCI 2009/2010 that I found in the CCPL
/// floor sum shenanigans
// taken from https://asfjwd.github.io/2020-04-24-floor-sum-ap/
long long FloorSumAP(long long a, long long b, long long c, long long n){}
// calculating $\sum_{x=0}^n \lfloor (ax+b)/c \rfloor$
  if(!a) return (b / c) * (n + 1);
  if(a >= c || b >= c) return ( ( n * (n + 1) ) / 2) * (a / c) + (n + 1) * (b / c) + FloorSumAP(a % c, b % c, c, n);
  long long m = (a * n + b) / c;
  return m * n - FloorSumAP(c, c - b - 1, a, m - 1);
}