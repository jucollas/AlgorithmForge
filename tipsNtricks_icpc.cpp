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

// ham-path
// phi(nd,msk)= strt\in N(nd) OR OR_{e\in N(nd)} phi(nd,msk|(1<<nd))
// In cases where I'm only looking for cycles containing msk, then I will start
// from every node; however, I don't need to erase the memory. The reason is that
// if such cycle exists I can find it with any strt\in msk, if there is none it won't
// matter which strt\in msk I choose.
// This idea can be seen in https://codeforces.com/problemset/problem/1804/E



/////// Range Xor Shenanigans (assuming range is [l,r])
// Question: ¿Is subset in range interesting? ¿Is it doable in O(1)?
// I used this in industrial nim and in https://codeforces.com/problemset/problem/2056/F2

// calculating xor in range
int xor_range(int m){
	//calculates XOR_{x=0}^m x in O(1). Observe it includes m.
	int res=(m&1)?0:m;
	return res^=((m+1)>>1)&1;
}
 
// The idea is that if I fix the bits that are always set, then I can exclude them
// x= 0010100
// y= ??1?1?? 
// meaning that if I exclude the 1's it becomes a xor_range problem
// this involves calculating the new m (since the actual m may not me superset)
//    and re-adding again the 'erased' bits after xor_range
int super_set(int m,int x){
	//calculates XOR_{y=0}^m (y&x==x)*y in O(1)
	if(x>=m)return 0;
	
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