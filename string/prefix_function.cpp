/* Autor: Oscar Vargas Pabon
Probado en mis biblioteca
*/ vector<int> prefix_function(const string&cad){//O(n)
// pi[i] -> el tamano del mayor prefijo 
// de cad que tambien es sujifo de cad[0..i]
	const int n = cad.size(); vector<int> pi( n, 0 );
	for(int i=1;i<n;++i){ pi[i] = pi[i-1];
		while(pi[i]&&cad[i]!=cad[pi[i]])pi[i]=pi[pi[i]-1];
		if(cad[i]==cad[pi[i]])++pi[i];
	} return pi;
} const int let='z'-'a'+1; // este fue probado en 
// https://codeforces.com/problemset/problem/1721/E
vector<vector<int>> prefix_automata(const string&cad){//$O(n*\Sigma)$
	vector<int> pi=prefix_function(cad);
	vector<vector<int>> fpi(cad.size(),vector<int>(let));
	for(int i=0;i<int(cad.size());++i)for(int j=0;j<let;++j){
        if(i>0&&cad[i]-'a'!=j)
             fpi[i][j]=fpi[pi[i-1]][j];
        else fpi[i][j]=i+(cad[i]-'a'==j);
    }return fpi; }