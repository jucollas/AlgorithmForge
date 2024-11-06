/*
Autor: Oscar Vargas Pabon
Lo probe en UVA 1449

La implementacion asume reinicializar 'states' - states=vector<Node>(1);
	* Notar que pattern no considera el vertice actual como un posible pattern
*/

const int K = 26;

class Node {
	public:
	int next[K];
	int output, parent, link, pattern;
	char letter;
	
	Node ( int p=0, char l=0, int o=-1 ) : parent(p), letter(l),output(o) {
		memset( next, -1, sizeof(next) );
		link=pattern=-1;
	}
};

vector<Node> states(1);

void add_string( const string & str, int strInd=0 ) {
	/* añade la cadena 'str' al trie, la marca con 'strInd' */
	int act = 0;
	for ( const char letter : str ) {
		char ch = letter-'a';
		if ( states[act].next[ch] == -1 ) {
			states[act].next[ch] = states.size();
			states.push_back( Node( act, letter ) );
		}
		act = states[act].next[ch];
	}
	states[act].output = max( states[act].output, strInd );
}

int link( int act ) ;

int next( int act, char letter ) {
	/* hallo el siguiente estado a visitar segun el actual y el caracter */
	char ch = letter-'a';
	if ( states[act].next[ch] == -1 ) {
		if ( act ) {
			states[act].next[ch] = next( link(act), letter );
		} else states[act].next[ch] = 0;
	}
	return states[act].next[ch];
}

int link( int act ) {
	/* Halla el indice al nodo que representa el 'longest proper suffix' en el trie */
	if ( states[act].link == -1 ) {
		if ( act && states[act].parent ) {
			states[act].link = next(link(states[act].parent),states[act].letter);
		} else states[act].link = 0;
	}
	return states[act].link;
}

int pattern( int act ) {
	/* Halla el siguiente elemento que pertenece a un patron */
	if ( states[act].pattern == -1 ) {
		int nxt = link( act );
		
		if ( !nxt ) states[act].pattern = 0;
		else if ( states[nxt].output != -1 ) states[act].pattern = nxt;
		else states[act].pattern = pattern( nxt );
	}
	return states[act].pattern;
}
