/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Probado en Codeforces -> ITMO -> SuffixArray -> step4 -> A
*/

const char EOS = '$';
// end-of-string -> se supone que es un caracter 
//           estrictamente menor a todos los demas del string

vector<int> suffixArray( string &str ) {
/* Halla el suffix-array (SA) en tiempo O(nlg n) */
    str.pb( EOS ); // ahorra edge-cases

    int i, n = str.size();
    vector<int> code( n ), master( n ), newCode( n );
    {
		// para las 'equivalence classes' cuando solo contienen un caracter
        vector<pair<char,int>> previousMaster( n ); 
        for ( i = 0 ; i < n ; ++i ) previousMaster[i] = pair<char,int> ( str[i], i );
        sort( previousMaster.begin(), previousMaster.end() );
		
        master[previousMaster[0].second] = 0; // ahorra edge-cases
        for ( i = 1 ; i < n ; ++i ) {
            if ( previousMaster[i-1].first < previousMaster[i].first )// distinto al anterior
                code[previousMaster[i].second] = code[previousMaster[i-1].second]+1;
            else // igual al anterior, mantiene el codigo
                code[previousMaster[i].second] = code[previousMaster[i-1].second];
        }
		// actualiza el master (SA computado hasta ahora)
        for ( i = 0 ; i < n ; ++i ) master[i] = previousMaster[i].second;
    }

    int k = 1;
    while ( k < n && code[master.back()] < n-1 ) {
		// hace el 'cyclic-shift' del master
        for ( i = 0 ; i < n ; ++i ) master[i] = (master[i]-k+n)%n;

        {
            vector<int> copy = master; // hago el bucket-sort en O(n)
            vector<int> bucket_size( n+1, 0 );
            for ( i = 0 ; i < n ; ++i ) ++bucket_size[code[copy[i]]+1];// cuento los elementos de cada 'bucket'
            for ( i = 1 ; i <= n ; ++i ) bucket_size[i] += bucket_size[i-1];//hago el arreglo de prefijos
            for ( i = 0 ; i < n ; ++i ) {
                master[bucket_size[code[copy[i]]]++] = copy[i];//lleno master otra vez
            }
        }

        newCode[master[0]] = 0; // creo las nuevas 'equivalence clases'
        for ( i = 1 ; i < n ; ++i ) {
            newCode[master[i]] = newCode[master[i-1]];
            if ( code[master[i-1]] != code[master[i]] || code[(master[i-1]+k)%n] != code[(master[i]+k)%n] ) {
                ++newCode[master[i]]; // la tupla de codigos es distinta
            }
        }
        code = newCode;
        k = ( k << 1 );
    }

    str.ppb(); // quito el EOS
    return master;
}

vector<int> lcpArray( vector<int> &sufix, string &str ) {
/* halla el longest-common-prefix array del sufixArray en tiempo O(n)
	Retorna lcp[i] -> lcp de sufix[i] y sufix[i+1] */
    str.push_back( EOS ); // me ahorra un edge-case
	int n = sufix.size();
	
    vector<int> inv ( n ); // para obtener la poscicion de i en el SA
    for ( int i = 0 ; i < n ; ++i ) inv[sufix[i]] = i;
	
    vector<int> lcp ( n-1 );
    int k = 0, j;
    for ( int i = 0 ; i < n-1 ; ++i ) {
        j = sufix[inv[i]-1]; // inv[i]-1 no se sale porque sufix[0] representa EOS
        while ( str[i + k] == str[j + k] ) ++k;
        lcp[inv[i]-1] = k;
        k = max( k-1, 0 );
    }
    str.pop_back(); // quito EOS
    return lcp;
}