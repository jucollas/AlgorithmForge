/* Tomado de https://cp-algorithms.com/string/manacher.html
Calcula p[i]=x donde [i-x,i+x] es un palindromo
Notar que no lo calculamos sobre la cadena, sino sobre
    c'='&'+a[0]+'#'+...+a[i]+...+'#'+a[n]+'&'
Si c'[i]='#' me refiero a un palindromo par (de lo contrario un impar) */
vector<int> manacher_odd(string s) {int n = s.size();
    s = "*" + s + "^"; vector<int> p(n + 2);
    int l = 0, r = 1; for(int i = 1; i <= n; i++) {
        if(i <= r) p[i] = min(r - i, p[l + (r - i)]);
        while(s[i - p[i]] == s[i + p[i]]) p[i]++;
        if(i + p[i] > r) l = i - p[i], r = i + p[i];
    } for(int i=0;i+1<p.size();++i)p[i]=p[i+1];
    p.pop_back();p.pop_back(); return p;
}vector<int> manacher(const string&s) {
    string t;for(auto c: s) t += string("|") + c;
    auto res = manacher_odd(t + "|");
    for(int i=0;i+1<res.size();++i)res[i]=res[i+1]-1;
    res.pop_back();if(!res.empty())res.pop_back();
    return res;}