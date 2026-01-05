// author: Oscar Vargas Pabon
// check for more info on this. The problem that made me do this;
// and it has some stolen impl with better times on that problem
// https://qoj.ac/contest/1814/problem/300

typedef complex<double> cmplx;
cmplx knroot(int n,int k){double tmp=double(2*k)*M_PI/double(n);return {cos(tmp),sin(tmp)};}
struct fft_info{
	vector<cmplx> omega,iomega;
	fft_info()=default;
	fft_info(int lgi){
		const int n=1<<lgi;
		omega.resize(n);iomega.resize(n);
		cmplx w=knroot(n,1),iw=knroot(n,n-1);
		omega[0]=iomega[0]=1;
		rep(i,1,n)omega[i]=omega[i-1]*w,iomega[i]=iomega[i-1]*iw;
	}
};array<fft_info,24> f_info;
void fft( vector<cmplx> &A,bool invert ){
	// Coley-Tuckey
	const int n=sz(A),lgi=ilog2(n);
	for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1; // counting in reverse
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(A[i], A[j]);
    } if(f_info[lgi].omega.empty())f_info[lgi]=fft_info(lgi);
	vector<cmplx> &om=(invert)?f_info[lgi].iomega:f_info[lgi].omega;
	for(int m=1,acb=1;m<=lgi;++m,acb<<=1){
		int l=1<<(lgi-m);
		for(int s=0;s<(1<<lgi);s+=1<<m){
			int ind=0; for(int j=0;j<acb;++j){
				cmplx v=A[s+j+acb]*om[ind];
				A[s+j+acb]=A[s+j]-v;
				A[s+j]+=v;
				
				ind=(ind+l)&(n-1);//scared of non-optimized mod
			}
		}
	}
}