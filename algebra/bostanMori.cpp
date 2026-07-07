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
mint familyFriendly_bostanMori( lint k, fps P, fps Q ){
	// computes |x^k|P/Q in time O(d*lgd*lgk)
	const int d=Q.size(); while(k){ 
		fps nQ=Q;rep(i,0,nQ.size())if(i&1)nQ[i]=-nQ[i];
		P*=nQ; Q*=nQ;
		
		rep(i,0,d)Q[i]=Q[i*2];
		Q.trunc(d);
	
		// P[d*2] to prevent issues of P[ind1]=P[ind2] when sz(P)>=ind2
		P[d*2]; rep(i,0,d)P[i]=P[i*2+(k&1ll)];
		P.trunc(d); k>>=1;
	} return P[0]/Q[0]; // the method below gives a rough ~1/4 speedup
} mint bostanMori( lint k, fps P, fps Q ){
    // computes |x^k|P/Q in time O(d*lgd*lgk)
    const int d=Q.size(),d2=1<<(ilog2(d-1)+2);
    vector<mint>&Pf=P.F,&Qf=Q.F,nQ(d2);
    Pf.resize(d2,0);Qf.resize(d2,0); while(k){ 
        rep(i,0,d2)nQ[i]=(i&1)?-Qf[i]:Qf[i];

        internal::fft(Qf,0);internal::fft(Pf,0);internal::fft(nQ,0);
        rep(i,0,d2)Pf[i]*=nQ[i],Qf[i]*=nQ[i];//saves 3 ffts
        internal::fft(Qf,1);internal::fft(Pf,1);
        
        rep(i,0,d)Q[i]=Q[i*2],P[i]=P[i*2+(k&1ll)];
        k>>=1; rep(i,d,d2)Q[i]=P[i]=0;
    } return P[0]/Q[0];
} mint seq_bostanMori(lint k,fps P,fps Q){
    // P are assumed to be the first |P| terms of the sequence
    // Q is assumed to be the characteristic polynomial of degree d
    // $F_i=\sum_j F_{i-j}q_j$ then $Q=x^d-\sum_{j=0}^{d-1} q_jx^j$
    //reverse(Q.begin(),Q.end());
    // the reverse transforms to representation
    // $Q=1-\sum_{j=0}^{d-1}x^{j+1}q_j$
    const int d=Q.size();
    P=(P.trunc(d)*Q).trunc(d-1);
    return bostanMori(k,P,Q); }