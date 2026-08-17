/* Author: Oscar Vargas Pabon
tested in https://codeforces.com/problemset/problem/1725/E

La mayor idea de esto es lo de que para todo $A\subseteq V$ se cumple que $H=\{lca(u,v)|u,v\in A\}$
	es igual a yo ordenar por tiempo e entrada del dfs y solo hacer $\{lca(A_i,A_{i+1})\}_{i=0}^{|A|-1}$

Asumo un lca y calculo de tiempos de entrada/salida
*/ // is u ancestor of v?
bool is_parent(int u,int v){return tin[u] < tin[v] && tout[u] >= tout[v];}
pair<vector<int>,vector<vector<int>>> build_virtualTree( const vector<int> &vertex) {
    //Retorna un arbol virtual en O(n lg n) con raiz en 0 donde n=|vertex|
    if ( vertex.empty() ) return make_pair(vector<int>(),vector<vector<int>>());
    const int n = vertex.size();
        
    vector<pair<int,int>> vtNode(n); // para hacer un sort por los tin de los vertices
    for(int i=0;i<n;++i)vtNode[i]=make_pair(tin[vertex[i]],vertex[i]);
    sort(vtNode.begin(),vtNode.end());

    for(int i=1;i<n;++i){ // añadir los LCA necesarios para especificar el arbol
        const int new_vertex=lca(vtNode[i-1].second,vtNode[i].second);
        vtNode.emplace_back(tin[new_vertex],new_vertex);
    } sort(vtNode.begin(),vtNode.end());

    vtNode.emplace_back( -1, -1 );// implementation vertex
        
    vector<vector<int>> tree(1,vector<int>());
    vector<int> ren{vtNode.front().second};
    tree.reserve(n);ren.reserve(n);
    
    vector<int> stack{0}; // el stack para construir el vt
    for(auto[tin_nd,tree_nm]:vtNode)if(tree_nm!=ren.back()){
        const int vt_nm=ren.size(); // there are cases where nodes get duplicated
        ren.push_back(tree_nm);
        tree.push_back( vector<int>() );
            
        while(int(stack.size())>=2&&(tree_nm==-1||!is_parent(ren[stack.back()],tree_nm))){
            const int tmp = stack.back(); stack.pop_back();
            // añado la arista entre este y el anterior elemento
            tree[tmp].push_back( stack.back() );
            tree[stack.back()].push_back( tmp );
        } stack.push_back(vt_nm); // para no contar el 'implementationVertex'
    } tree.pop_back();ren.pop_back(); return {ren,tree};
}