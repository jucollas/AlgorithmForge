/*
Autor: Oscar Vargas Pabon

Lo desarrolle para: 11476 Factorizing Large Integers
Notar que nunca obtuve AC por este metodo (es muy lento pero interesante)

idea tomada de:
	"Acceleration of Wheel Factoring Techniques"
	Alaa M. Zaki 1 , M. E. Bakr 2, Arwa M. Alsahangiti 2, Saima Khan Khosa 3 and Khaled A. Fathy 4,*
	
Invocar 'precompute' en el main

precompute(1e4, template_limit);// 8 primos -> 17.1024% del espacio de busqueda
*/

const int template_limit = 1e7;
int cribe[template_limit], primes[template_limit];
int wheel[template_limit], wheel_sz, prime_sz, init;

void precompute( int n, int p_max ) {
	memset( cribe, -1, sizeof(cribe) );
	cribe[0] = cribe[1] = 1;
	int ind = 2;
	while ( ind < n ) {
		cribe[ind] = ind;//criba de toda la vida (version extrana mia)
		for ( int i = ind*ind ; i < n ; i+= ind ) cribe[i] = ind;
		while ( ind < n && cribe[ind]>0 ) ++ind;
	}
	
	int val = 1; prime_sz = 0; // sacar los primos y el modulo de la wheel
	for ( int i = 2 ; val*i < p_max && i < n ; ++i ) {
		if ( cribe[i] == i ) {
			primes[prime_sz++] = i;
			val *= i;
		}
	}
	init = primes[prime_sz-1]+2;// valor inicial para la wheel
	
	int lst = init; wheel_sz=0;
	for ( int i = 0 ; i < val ; ++i ) {
		bool plaus = 1; // reviso si pertenece a la rueda
		for ( int j = 0 ; j < prime_sz && plaus ; ++j ) plaus = (init+i)%primes[j]>0;
		if ( plaus ) { // lo agrego a la rueda
			wheel[wheel_sz++] = i+init-lst;
			lst=i+init;
		}
	}
	wheel[wheel_sz++] = init+val+wheel[0]-lst;//cierro el ciclo
	if ( !wheel[0] ) {
		// edge-case
		for ( int i = 1 ; i < wheel_sz ; ++i ) wheel[i-1] = wheel[i];
		--wheel_sz;
	}
}


map<int,int> wheel_fact(int x) {
	map<int,int> decomp;
	for ( int i = 0 ; i < prime_sz ; ++i ) {
		int p = primes[i];
		while ( n%p==0ll ) {
			n /= p; ++decomp[p];
		}
	}

	int p = init, ind = 0;
	while ( n>1ll && p*p <= n ) {
		while ( n%p==0ll ) {
			n/=p; ++decomp[p];
		}
		
		p += wheel[ind++];
		if ( ind == wheel_sz ) ind = 0;
	}
	if ( n>1ll ) ++decomp[n];// n es primo
	return decomp;
}
