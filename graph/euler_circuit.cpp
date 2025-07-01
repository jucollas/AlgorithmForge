/*
Autor: Oscar Vargas Pabon

Esta implementacion destruye el grafo y es solo para dirigidos.

Recordar que solo hay un circuito euleriano si indeg[nd]==outdeg[nd] para
	todos los nodos.
Un camino euleriano se cumple cuando para u,v indeg[u]+1==outdeg[u],
	indeg[v]==outdeg[v]+1 y todos los demas nodos cumplen indeg[nd]==outdeg[nd].
	En este caso basta con iniciar la primera solucion de 'euler_circuit' en u.
*/
void euc_aux( vector<list<pair<int,int>>> &g, int nd, list<pair<int,int>> &res ) {
	if ( !g[nd].empty() ) {
		int neo_nd,val;
		neo_nd = g[nd].front().first; val = g[nd].front().second;
		
		res.pb( pair<int,int>(nd,val) );
		g[nd].ppf();
		
		euc_aux( g, neo_nd, res );
		
	}
}

list<pair<int,int>> euler_circuit( vector<list<pair<int,int>>> &g ) {
	list<pair<int,int>> res;
	euc_aux( g, 0, res );
	for ( list<pair<int,int>>::iterator it = res.begin() ; it != res.end() ; ++it ){
		if ( !g[it->first].empty()  ) {
			list<pair<int,int>> tmp;
			euc_aux( g, it->first, tmp );
			
			for ( const pair<int,int> &act : tmp ) res.insert(it,act);
			while ( !tmp.empty() ) {
				tmp.ppf(); --it;
			}
		}
	}
	return res;
}

