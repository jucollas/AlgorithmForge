template<typename tint> constexpr int countr_zero_constexpr(tint n){if(n==tint(0))return -1;int x=0;while(!(n&(tint(1)<<x)))++x;return x;}
template<typename tint> constexpr tint primitive_root(tint m){
    if ((27ull<<59)+1 == m ) return 5;
    if( (549755813881ull<<24)+1 == m || 998244353 == m )return 3;
    tint dec[32]={},dind=0,md=m-1; for(tint i=2;i*i<=md;++i)if(md%i==0){
        dec[dind++]=i; while(md%i==0)md/=i;
    } if(md>1)dec[dind++]=md; //tests for g such that $\forall_{p|(md-1)} g^{(md-1)/p}=1(mod md)$
    bool fnd=0; tint pr=1;while(!fnd){ ++pr; fnd=1; 
        for(tint i=0;i<dind&&fnd;++i){ // by properties, this will always end
            if constexpr ( std::is_same_v<tint,unsigned long long> )
                fnd=mpow<tint,__uint128_t>(pr,(m-1)/dec[i],m)!=1;
            else if constexpr ( std::is_same_v<tint,long long> )
                fnd=mpow<tint,__uint128_t>(pr,(m-1)/dec[i],m)!=1;
            else fnd=mpow<tint>(pr,(m-1)/dec[i],m)!=1; }
    } return pr;// std::cerr <<pr << " _ primitive root" << endl;
} template<typename tfps> constexpr auto prec_rank_fft_root() {
    constexpr int rank=countr_zero_constexpr(tfps::mod()-1);
    std::array<tfps,rank+1> root; // precompute the roots
    root[rank]=tfps(primitive_root(tfps::mod())).pow((tfps::mod()-1)>>rank);
    for(int i=rank;i>0;--i)root[i-1]=root[i]*root[i];
    return root;
}template<typename tfps> inline constexpr auto rank_fft_root = prec_rank_fft_root<tfps>();
template<typename tfps,uint64_t N> inline constexpr tfps fft_root(){
    if constexpr(N<=0ll||countr_zero_constexpr(N)>=int(rank_fft_root<tfps>.size()))return tfps(0);
    return rank_fft_root<tfps>[countr_zero_constexpr(N)];
} template<typename tfps,int n>void butterfly(tfps*a){
    constexpr int m=n/2; if constexpr(n<=1)return;
    butterfly<tfps,m>(a);butterfly<tfps,m>(a+m);
    constexpr tfps wlen=fft_root<tfps,n>();tfps w=1;
    for(int i=0;i<m;++i){tfps u=a[i],v=w*a[i+m];
        a[i]=u+v; a[i+m]=u-v; w*=wlen; }
}template<typename tfps,int n=1<<25>void fft(vector<tfps>&a,bool sd){
    if constexpr(n<=1)return;
    if(n>a.size())fft<tfps,n/2>(a,sd);
    else{ for(int i=1,j=0;i<n;++i){ int bit=n>>1;
            for(;j&bit;bit>>=1)j^=bit;//bit reversal
            j^=bit; if(i<j)swap(a[i],a[j]);
        } butterfly<tfps,n>(a.data()); if(sd){
            tfps ni=tfps(n).inv();//multiplicative
            reverse(a.begin()+1,a.end());// inverse
            for(tfps&ac:a)ac*=ni; } }
} template<typename tfps> vector<tfps>convolution(vector<tfps>A,vector<tfps>B){
    int n=A.size()+B.size()-1,m=1;while(m<n)m*=2;
    A.resize(m,0);B.resize(m,0);
    fft<tfps>(A,0);fft<tfps>(B,0);
    for(int i=0;i<m;++i)A[i]*=B[i];
    fft<tfps>(A,1);A.resize(n);return A;
} template<typename tfps>vector<tfps>inverse(vector<tfps>F,int n){
    assert(F[0]);// G*F=1
    vector<tfps>G={F[0].inv()}; for(int e=2;e<2*n;e<<=1){
        vector<tfps> ac=convolution(G,// gives a ~/2 speedup
                {F.begin(),F.begin()+min<int>(F.size(),e)}
        ); for(tfps&act:ac)act=-act;
        ac[0]+=2;G=convolution(G,ac);
        G.resize(e);
    } G.resize(n);return G; }