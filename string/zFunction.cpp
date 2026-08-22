/* Autor: Oscar Vargas Pabon
Lo probe en mi carpeta de pruebas
*/ vector<int> zFunction(const string&cad){ // O(n). 
//  z[i] -> maximo prefijo comun de cad y cad[i...] 
	const int n = cad.size(); vector<int> z(n,0);
	int l=-1,r=-1;for(int i=1;i<n;++i){ z[i]=r-i;
		if(z[i]>z[i-l])z[i]=z[i-l];if(z[i]<0)z[i]=0;
		while(i+z[i]<n&&cad[z[i]]==cad[i+z[i]])++z[i];
		if(i+z[i]>r)r=i+z[i],l=i;
	} return z; }