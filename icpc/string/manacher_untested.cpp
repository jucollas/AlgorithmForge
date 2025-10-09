/*
Autor: Oscar Vargas Pabon
Tomado de https://cp-algorithms.com/string/manacher.html
No ha sido testeado

Calcula p[i]=x donde [i-x,i+x] es un palindromo
Notar que no lo calculamos sobre la cadena, sino sobre
	c'='&'+a[0]+'$'+...+a[i]+...+'$'+a[n]+'&'
Si c'[i]='$' me refiero a un palindromo par (de lo contrario un impar)
*/

vector<int> manacher( const string &cad ){
	string c2="&$";
	for(char c:cad){c2.push_back(c);c2.push_back('$');}
	c2.push_back('%');
	int n=c2.size();
	vector<int> res(n-1,-1);
	int l=0,r=1;
	for(int i = 1 ; i < n-1 ; ++i ) {
		res[i]=min(r-i,res[l+(r-i)]);
		while(c2[i-res[i]]==c2[i+res[i]])++res[i];
		if(i+res[i]>r)l=i-res[i],r=i+res[i];
	}
	return res;
}