/* Author: Oscar vargas Pabon
Tested in: https://www.luogu.com.cn/problem/P4929
Guided by: https://www.luogu.com.cn/article/uz2ommgq
            https://www.luogu.com.cn/article/drdgqy2s
            https://www.cnblogs.com/JMXZ/articles/17816659.html
A really weird rabbit hole to go into. The order of the walk when erasing and when
    inserting must be inverted, as it may happen that otherwise, certain nodes
    end up with bad values in them. So they are important.
It is also widely assumed in my impl that the 0 is the index of the 'row/col' intersection.
    Furthermore, the rows indexes are 1..n; col indexes are n+1..m+n+1
There are other applications of the idea of DLX, yet I wont go into them as of now.
*/
const int max_n=505,max_o=5005+2*max_n+33; struct DLX{
    int u[max_o],d[max_o],r[max_o],l[max_o],n,m,am;//for nodes
    inline void ins(int x){ // inserts node x
        u[d[x]]=d[u[x]]=l[r[x]]=r[l[x]]=x;
    } void init(int nn,int mm){
        // special-node index |-> 0
        // special row-index |-> 1..n
        // special col-index |-> n+1..m+n+1
        n=nn;m=mm;am=n+m+1;u[0]=d[0]=r[0]=l[0]=0;
        rep(i,1,n+1){ l[i]=r[i]=i;u[i]=i-1;d[i]=0; ins(i); }
        rep(i,n+1,n+m+1){u[i]=d[i]=i;l[i]=i-1==n?0:i-1;r[i]=0; ins(i); }
    } int ares[max_n],ram,cl_am[max_n],cl_nm[max_o];
    bool solve(){ memset(cl_am,0,sizeof(int)*(m+1));
        {int cl=0;do{ // cl_am and cl_nm filled
            for(int rw=d[cl];rw!=cl;rw=d[rw]) ++cl_am[cl_nm[rw]=cl?cl-n:0];
            cl=r[cl]; } while(cl); }
        ram=0; auto aux=[&](auto rec)->bool{ if(!r[0]||!cl_am[0])return !r[0];
            int piv=-1;for(int cl=r[0];cl;cl=r[cl]){
                if(!cl_am[cl-n])return 0;//check column-amount heuristic
                if(piv==-1||cl_am[cl-n]<cl_am[piv-n])piv=cl;//for some 
            } for(int rw=d[piv];rw!=piv;rw=d[rw]){// reason taking the min
                for(int cl=r[rw];cl!=rw;cl=r[cl]){//is better
                    u[d[cl]]=u[cl];d[u[cl]]=d[cl];
                    --cl_am[cl_nm[cl]];
                }//remove pivot column and associates
            } l[r[piv]]=l[piv];r[l[piv]]=r[piv];
            bool res=0; for(int nd=d[piv];nd!=piv&&!res;nd=d[nd]){
                for(int rw=l[nd];rw!=nd;rw=l[rw]){//the ones in same row as nd
                    if(rw<=n){ ares[ram++]=rw;
                        u[d[rw]]=u[rw];d[u[rw]]=d[rw];
                        continue;
                    } const int mrw=cl_nm[rw]+n; l[r[mrw]]=l[mrw],r[l[mrw]]=r[mrw];
                    for(int cl=d[mrw];cl!=mrw;cl=d[cl])// the ones in same col as rw
                        for(int x=r[cl];x!=cl;x=r[x])u[d[x]]=u[x],d[u[x]]=d[x];
                } res=rec(rec); ram-=!res; if(res)return 1;
                for(int rw=r[nd];rw!=nd;rw=r[rw]){//the ones in same row as nd
                    if(rw<=n){ d[u[rw]]=u[d[rw]]=rw; continue;}
                    const int mrw=cl_nm[rw]+n; l[r[mrw]]=r[l[mrw]]=mrw;
                    for(int cl=u[mrw];cl!=mrw;cl=u[cl])// the ones in same col as rw
                        for(int x=l[cl];x!=cl;x=l[x]) u[d[x]]=d[u[x]]=x;
                }
            } for(int rw=u[piv];rw!=piv;rw=u[rw]){
                for(int cl=l[rw];cl!=rw;cl=l[cl]){
                    u[d[cl]]=d[u[cl]]=cl; ++cl_am[cl_nm[cl]];
                }//insert pivot column and associates
            } l[r[piv]]=r[l[piv]]=piv; return res;
        }; return aux(aux); }
} dlx;
void solve() {
    int n,m;cin>>n>>m;
    dlx.init(n,m); rep(i,0,n){
        int lind=i+1;rep(j,0,m){
            int x;cin>>x;if(!x)continue;
            const int aind=dlx.am++;
            dlx.l[aind]=lind;dlx.r[aind]=i+1;
            dlx.u[aind]=dlx.u[j+n+1]; dlx.d[aind]=j+n+1;
            dlx.ins(lind=aind);
        }
    } if(dlx.solve()){
        rep(i,0,dlx.ram)cout << dlx.ares[i] << " \n"[i+1==dlx.ram];
    } else cout << "No solution!\n";
}