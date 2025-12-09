/*
Author:Oscar Vargas Pabon

Taken from:
	* https://codeforces.com/blog/entry/111862
	* "A Simple and Fast Algorithm for Computing the N-th Term of a Linearly Recurrent Sequence"
		Alin Bostan and Ryuhei Mori
Note that other versions as the one computing a 'window' of terms are still not implemented

I assume fps from my impl

Notes: If I want to generate the first n terms then I can also do it by doing
			(P*Q.inv(n)).trunc(n); for F=P/Q the polinomials
			
*/

mint bostanMori( lint k, fps P, fps Q ){
	// computes kth element of the recurrence 
	// F=P/Q where Q is the characteristic polynomial
	// F_i=\sum_j F_{i-j}d_j then Q=x^{sz(d)}-\sum_j d_j
	
	const int d=sz(Q);reverse(all(Q));
	P=(P.trunc(sz(Q))*Q).trunc(sz(Q));
	
	while(k){
		fps nQ=Q;rep(i,0,sz(nQ))if(i&1)nQ[i]=-nQ[i];
		P*=nQ; Q*=nQ;
		
		rep(i,0,d)Q[i]=Q[i*2];
		Q.trunc(d);
	
		// P[d*2] to prevent issues of P[ind1]=P[ind2] when sz(P)>=ind2
		P[d*2]; rep(i,0,d)P[i]=P[i*2+(k&1ll)];
		P.trunc(d);
		
		k>>=1;
	}
	return P[0]/Q[0];
}