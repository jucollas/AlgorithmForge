/*
Autor: Oscar Vargas Pabon
These are various tricks to be used at will
*/

// iterating over ranges [l_i,r_i]\subseteq[L,R] such that (assuming integer division)
//////   \forall_{j,h\in[l_i,r_i]}n/j=n/h
/// can be seen in https://codeforces.com/problemset/problem/2072/G
void diviter(int n, int L, int R ){
	int l = L,r;
	while( l <= R ){
		int nh = n/l;
		r = min(R,n/nh);
		
		/// do stuff with the range i\in[l,r] n/i=nh
		
		l=r+1;
	}
}