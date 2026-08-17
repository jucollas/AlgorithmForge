/* Author: Oscar Vargas Pabon
Partially tested in https://qoj.ac/contest/1814/problem/843
It may be slightly slow
*/
typedef vector<vector<mint>> mat;
template<typename tint>mat transpose(const mat&m){
    const int n=m.size(),k=m[0].size();
    vector<vector<tint>>tran(k,vector<tint>(n));
    for(int i=0;i<n;++i)for(int j=0;j<k;++j)tran[j][i]=m[i][j];
    return tran;
} template<typename tint>mat mat_mult(const mat&m1,const mat&m2){
    const int n=m1.size(),m=m2.size(),k=m2[0].size();
    assert(int(m1[0].size())==m);
    vector<vector<mint>> res(n,vector<mint>(k,tint(0)));
    for(int i=0;i<n;++i)for(int j=0;j<m;++j)for(int h=0;h<k;++h)res[i][h]+=m1[i][j]*m2[j][h];
    return res;
} template<typename tint>tint trace(const mat&m){
    tint res=0;for(int i=0;i<int(m.size());++i)res+=m[i][i];
    return res;
}template<typename tint>tint determinant(mat m){
    const int n=m.size();assert(int(m[0].size())==n);
    tint det=1;for(int i=0;i<n;++i){
        for(int j=i;j<n;++j)if(m[i][j]!=tint(0)){
            swap(m[i],m[j]); if(j!=i)det=-det;
            break;
        }det*=m[i][i];if(det==tint(0))return tint(0);
        {tint ac=m[i][i].inv();for(int j=i;j<n;++j)m[i][j]*=ac;}
        for(int j=i+1;j<n;++j)if(m[j][i]!=tint(0)){
            tint ac=-m[j][i];for(int k=i;k<n;++k)m[j][k]+=m[i][k]*ac;
        }
    }return det;
}template<typename tint>
tint kth_lin_recurrence(uint64_t e,const vector<tint>&recur,const vector<tint>&init){
// I assume the recurrence coefficients to come as
// $a(n)=\sum_{i=0}^{n-1}recur(n-1-i)a(n-1-i)$ for the sequence $a$
    const int n=recur.size();assert(int(init.size())==n);
    if(e<n)return init[e];
    vector<vector<tint>> A(n,vector<tint>(n,tint(0))),B=A;
    for(int i=0;i<n;++i)B[i][i]=1;
    for(int i=1;i<n;++i)A[i-1][i]=1;
    A.back()=recur; e-=n; while(e){
        if(e&1ull)B=mat_mult<tint>(B,A);
        e>>=1;A=mat_mult<tint>(A,A);
    } return mat_mult<tint>(B,transpose<mint>(vector<vector<tint>>(1,init))).back().back();
}