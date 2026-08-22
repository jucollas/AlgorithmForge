/* Autor: Oscar Vargas Pabon
Probado en Codeforces -> ITMO -> SuffixArray -> step4 -> A
*/ const char EOS = '#';//end-of-string -> estrictamente menor
//                            a todos los demas del string
inline int m_add(int x,int y,int m){return m-x>y?x+y:y-(m-x);}
inline int m_sub(int x,int y,int m){return x>=y?x-y:x+(m-y);}
vector<int> suffixArray( string &str ) {
/* Halla el suffix-array (SA) en tiempo O(nlg n) */
    str.push_back( EOS ); const int n=str.size(); // ahorra edge-cases
    vector<int> master(n),code(n),newCode(n); {
        // para las 'equivalence classes' cuando solo contienen un caracter
        vector<pair<char,int>> p_master(n); 
        for(int i=0;i<n;++i)p_master[i]=pair<char,int>(str[i],i);
        std::sort(p_master.begin(),p_master.end());
        master[p_master[0].second]=0;for(int i=1;i<n;++i){
            if(p_master[i-1].first<p_master[i].first)
                 code[p_master[i].second]=code[p_master[i-1].second]+1;
            else code[p_master[i].second]=code[p_master[i-1].second]  ;
        } for(int i=0;i<n;++i)master[i]=p_master[i].second;
    }vector<int> copy,bucket(n);
    int k=1;while(k<n&&code[master.back()]<n-1){ // hace el 'cyclic-shift'
        for(int i=0;i<n;++i)master[i]=m_sub(master[i],k,n); {//  del master
            copy=master;fill(bucket.begin(),bucket.end(),0);// bucket-sort
            for(int i=0;i<n;++i)bucket[code[i]+1]++;//<- equivalent to 
            for(int i=1;i<n;++i)bucket[i]+=bucket[i-1];// code[copy[i]+1]
            for(int i=0;i<n;++i)master[bucket[code[copy[i]]]++]=copy[i];
        } newCode[master[0]]=0; for(int i=1;i<n;++i){// creo las nuevas 
            newCode[master[i]] = newCode[master[i-1]];//equivalence clases
            if(code[      master[i-1]     ]!=code[      master[i]     ]|| 
               code[m_add(master[i-1],k,n)]!=code[m_add(master[i],k,n)] )
                ++newCode[master[i]]; // la tupla de codigos es distinta
        } code = newCode; k = ( k << 1 );
    } str.pop_back(); return master; // quito el EOS
} vector<int> lcpArray( vector<int> &sufix, string &str ) {
/* halla el longest-common-prefix array del sufixArray en tiempo O(n)
    Retorna lcp[i] -> lcp de sufix[i] y sufix[i+1] */
    // me ahorra un edge-case
    str.push_back( EOS ); const int n = sufix.size();
    vector<int> lcp(n-1),inv(n); for(int i=0;i<n;++i)inv[sufix[i]]=i;
    int k=0; for(int i=0;i<n-1;++i) {
        const int j=sufix[inv[i]-1];// inv[i]-1 no se sale porque 
        while(str[i+k]==str[j+k])++k;// sufix[0] representa EOS
        lcp[inv[i]-1] = k;
        --k;if(k<0)k=0;
    } str.pop_back(); return lcp; } // quito EOS