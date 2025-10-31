/*
Autor: Oscar Vargas Pabon

GCD_CONV Currently UNTESTED
LCM_CONV tested in https://codeforces.com/gym/105053/problem/G

Note it wont work well for A*A

Taken from https://codeforces.com/blog/entry/119082

Im assuming from my template:
#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
*/


template<typename T>
class div_conv{
	vector<int> prm;
public:
	div_conv(int limit){
		vector<bool> crb(limit,1);crb[0]=crb[1]=0;
		rep(i,0,limit)if(crb[i]){
			prm.push_back(i);
			if(i*1ll*i<lint(limit))for(int j=i*i;j<limit;j+=i){
				crb[j]=0;
			}
		}
	}

	void zeta_div(vector<T> &vc){
		// vc'_x=\sum_{y|x}vc_y
		for(int p:prm)for(int i=1;i*p<int(vc.size());++i){
			vc[i*p]+=vc[i];
		}
	}
	void mobi_div(vector<T> &vc){
		// vc_x=\sum_{y|x}vc'_y
		for(int p:prm)for(int i=(vc.size()-1)/p;i;--i){
			vc[i*p]-=vc[i];
		}
	}
	void gcd_conv(vector<T> &A,vector<T> &B){
		// I assume |A|=|B|
		// A'_x=\sum_{x=gcd(u,v)}A_uB_v
		zeta_div(A);zeta_div(B);
		rep(i,0,A.size())A[i]*=B[i];
		movi_div(A);
	}
	
	void zeta_mult(vector<T> &vc){
		// vc_x=\sum_{x|y}vc'_y
		for(int p:prm)for(int i=1;i*p<int(vc.size());++i){
			vc[i*p]+=vc[i];
		}
	}
	void mobi_mult(vector<T> &vc){
		// vc'_x=\sum_{x|y}vc_y
		for(int p:prm)for(int i=(vc.size()-1)/p;i;--i){
			vc[i*p]-=vc[i];
		}
	}
	
	void lcm_conv(vector<T> &A,vector<T> &B){
		// I assume |A|=|B|
		// A'_x=\sum_{x=lcm(u,v)}A_uB_v
		zeta_mult(A);zeta_mult(B);
		rep(i,0,A.size())A[i]*=B[i];
		mobi_mult(A);
	}

};

div_conv<double> dc(template_limit);