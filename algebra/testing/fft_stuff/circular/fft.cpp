namespace internal{
constexpr int countr_zero_constexpr(unsigned int n) { int x = 0; while (!(n & (1 << x))) x++; return x; }
constexpr int primitive_root(int m){
	//tests for g such that $\forall_{p|(md-1)} g^{(md-1)/p}=1(mod md)$
	int dec[32]={},dind=0, md=m-1;
	for(int i=2;i*i<=md;++i)if(md%i==0){
		dec[dind++]=i; while(md%i==0)md/=i;
	} if(md>1)dec[dind++]=md;
	
	bool fnd=0; int pr=1;while(!fnd){
		++pr; fnd=1; // by properties, this will always end
		for(int i=0;i<dind&&fnd;++i)fnd=mpow<int>(pr,(m-1)/dec[i],m)!=1;
	} return pr;// std::cerr <<pr << " _ primitive root" << endl;
} template<typename tfps,int rank> constexpr tfps fft_root(){
    constexpr tfps c_root=primitive_root(tfps::mod());
    constexpr int rank2=countr_zero_constexpr(tfps::mod()-1);
    if constexpr(rank>=rank2)
        return c_root.pow((tfps::mod()-1)>>rank);
    else return fft_root<tfps,rank+1>()*fft_root<tfps,rank+1>();
} template<typename tfps,int rank> constexpr tfps fft_iroot(){
    constexpr int rank2=countr_zero_constexpr(tfps::mod()-1);
    if constexpr(rank>=rank2)
        return fft_root<tfps,rank>().inv();
    else return fft_iroot<tfps,rank+1>()*fft_iroot<tfps,rank+1>();
} template<typename tfps,int n> void butterfly_rec(tfps*a){
    if constexpr(n<=1)return;
    constexpr int e=ilog2(n),m=n/2;
    constexpr tfps wlen=fft_root<tfps,e>();
    tfps w=1; for(int i=0;i<m;++i){
        const tfps u=a[i],v=a[i+m];
        a[i  ]=u+v; a[i+m]=(u-v)*w;
        w*=wlen;
    }butterfly_rec<tfps,m>(a);butterfly_rec<tfps,m>(a+m);
}template<typename tfps,int n=(1<<countr_zero_constexpr(tfps::mod()-1))>
void butterfly(std::vector<tfps>&a){
    if constexpr(n<=1)return;
    if (n==int(a.size()))butterfly_rec<tfps,n>(a.data());
    else butterfly<tfps,n/2>(a);
} template<typename tfps,int n> void ibutterfly_rec(tfps*a){
    if constexpr(n<=1)return;
    constexpr int e=ilog2(n),m=n/2;
    constexpr tfps wlen=fft_iroot<tfps,e>();
    ibutterfly_rec<tfps,m>(a); ibutterfly_rec<tfps,m>(a+m);
    tfps w=1;for(int i=0;i<m;++i){
        const tfps u=a[i],v=w*a[i+m];
        a[i  ]=u+v; a[i+m]=u-v;
        w*=wlen; }
}template<typename tfps,int n=(1<<countr_zero_constexpr(tfps::mod()-1))>
void ibutterfly(std::vector<tfps>&a){
    if constexpr(n<=1)return;
    if (n==int(a.size()))ibutterfly_rec<tfps,n>(a.data());
    else ibutterfly<tfps,n/2>(a);
} template<typename tfps> void fft(vector<tfps>&a,bool invert){
    if(invert){ internal::ibutterfly(a);
        tfps n_1 = tfps(int(a.size())).inv();
        for (tfps&x:a)x*=n_1;
	}else internal::butterfly(a);
} template<typename tfps> void transposed_fft(std::vector<tfps> &A,bool invert){
	if (!invert) internal::ibutterfly<tfps>(A);
	reverse(A.begin() + 1, A.end()); if(invert){
		internal::butterfly<tfps>(A);
		tfps n_1=tfps(A.size()).inv();
		for (tfps &x: A) x *= n_1; }
}
}