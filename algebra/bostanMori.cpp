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
template<typename tfps> tfps bostanMori(std::vector<tfps>P,std::vector<tfps>Q,__uint64_t k){
    const int d=Q.size(),dlgi=ilog2(d-1)+2,d2=1<<dlgi,d22=d2/2;
    P.resize(d2,tfps(0));Q.resize(d2,tfps(0)); // computes |x^k|P/Q in time O(d*lgd*lgk)
    
    vector<tfps>wrot(d22);{ const tfps wlen=internal::rank_fft_root<tfps>[dlgi].inv();
        wrot[0]=1;for(int i=1;i<d22;++i)wrot[i]=wrot[i-1]*wlen;
        for(int i=0;i<d22;++i){
            int x=0;for(int e=0;e<dlgi-1;++e)x|=(i>>e&1)<<(dlgi-2-e);
            if(x<i)swap(wrot[i],wrot[x]);
        }for(int i=0;i<d22;++i)wrot[i]*=tfps(2).inv();
    } internal::fft(Q,0);internal::fft(P,0); while(k){
        for(int i=0;i<d2;++i)P[i]*=Q[i^1];  
        if(k&1ull) for(int i=0;i<d22;++i)P[i]=(P[i*2]-P[i*2+1])*wrot[i];
        else for(int i=0;i<d22;++i)P[i]=(P[i*2]+P[i*2+1])*tfps(2).inv();
        internal::fft_doubling(P.data(),d2);

        for(int i=0;i<d22;++i)Q[i]=Q[i<<1]*Q[i<<1|1];
        internal::fft_doubling(Q.data(),d2);
        k>>=1;
    } internal::fft(Q,1);internal::fft(P,1);
    return P[0]/Q[0];
} mint seq_bostanMori( fps P,fps Q,lint k){
    // P are assumed to be the first |P| terms of the sequence
    // Q is assumed to be the characteristic polynomial of degree d
    // $F_i=\sum_j F_{i-j}q_j$ then $Q=x^d-\sum_{j=0}^{d-1} q_jx^j$
    //reverse(Q.begin(),Q.end());
    // the reverse transforms to representation
    // $Q=1-\sum_{j=0}^{d-1}x^{j+1}q_j$
    const int d=Q.size();
    P=(P.trunc(d)*Q).trunc(d-1);
    return bostanMori<mint>(P.F,Q.F,k); }