/*
 * Hungarian Algorithm 
 * ----------------------------------------
 *  Entrada: 
 *   - Una matriz cuadrada cost_matrix de tamano nXn,
 *     donde cost_matrix[i][j] representa el costo de asignar
 *     el trabajador i al trabajo j.
 *
 *  Salida:
 *   - El costo minimo total de asignacion optima.
 *
 *  Uso:
 *   - Resuelve el "Assignment Problem": asignar n tareas a n agentes
 *     minimizando el costo total (o maximizando utilidad si se invierte el signo).
 *   - Complejidad: O(n^3)
 */

#include <bits/stdc++.h>
using namespace std;

class Hungarian {
private:
    vector<vector<int>> cost_matrix; // Matriz de costos
    int n;                           // Tamano de la matriz
    vector<int> u, v, p, way;        // Vectores auxiliares para el algoritmo

public:
    // Constructor: inicializa con la matriz de costos
    Hungarian(const vector<vector<int>> &matrix) : cost_matrix(matrix) {
        n = matrix.size();
        u.resize(n + 1, 0);
        v.resize(n + 1, 0);
        p.resize(n + 1, 0);
        way.resize(n + 1, 0);
    }

    // Funcion principal: devuelve el costo minimo total
    int compute() {
        for (int i = 1; i <= n; ++i) {
            p[0] = i;
            vector<int> minv(n + 1, INT_MAX);
            vector<bool> used(n + 1, false);
            int j0 = 0;
            while (true) {
                used[j0] = true;
                int i0 = p[j0], delta = INT_MAX, j1 = 0;
                for (int j = 1; j <= n; ++j) {
                    if (!used[j]) {
                        int cur = cost_matrix[i0 - 1][j - 1] - u[i0] - v[j];
                        if (cur < minv[j]) {
                            minv[j] = cur;
                            way[j] = j0;
                        }
                        if (minv[j] < delta) {
                            delta = minv[j];
                            j1 = j;
                        }
                    }
                }
                for (int j = 0; j <= n; ++j) {
                    if (used[j]) {
                        u[p[j]] += delta;
                        v[j] -= delta;
                    } else {
                        minv[j] -= delta;
                    }
                }
                j0 = j1;
                if (p[j0] == 0)
                    break;
            }
            while (true) {
                int j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
                if (j0 == 0)
                    break;
            }
        }

        int total_cost = 0;
        for (int j = 1; j <= n; ++j)
            total_cost += cost_matrix[p[j] - 1][j - 1];
        return total_cost;
    }
};
