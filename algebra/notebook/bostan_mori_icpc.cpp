// notice that some speedups like optimizing 1 fft of nQ
// or using recurrences among the even/odd terms of the FFT(F) may give
// greater speedups
template<typename tfps>tfps bostanMori(uint64_t k,vector<tfps>P,vector<tfps>Q){
    // computes [x^k]P/Q in time O(d*lgd*lgk)
    const int d=Q.size(); while(k){ vector<tfps> nQ=Q;
        for(int i=0;i<nQ.size();++i)if(i&1)nQ[i]=-nQ[i];
        P=convolution(P,nQ);Q=convolution(Q,nQ);
        for(int i=0;i<d;++i)Q[i]=Q[i*2];
        Q.resize(d,0); P.resize(d*2+1,0);
        for(int i=0;i<d;++i)P[i]=P[i*2+(k&1ll)];
        P.resize(d); k>>=1;
    } return P[0]/Q[0]; }