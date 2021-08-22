#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include "../Graphs/graph.h"

// Стуктура путь
struct Path{
    int l; // общая длина пути
    int *A; // массив вершин пути
};

void delete_path(struct Path *P){
    free(P -> A);
    free(P);
}

void output_path(struct Path *P, int V){
    printf("\nPath: %i ", P -> A[0]);
    for (int i = 0; i < V; i++)
        printf("--> %i ", P -> A[(i + 1) % V]);
    printf("\nPath length: %i\n", P -> l);
}

double probability(int i, int j, int *tabu, int V, int alpha, int beta, double **tau, double **eta){
    if (tabu[j])
        return 0.0;
    double P = pow(tau[i][j], alpha) * pow(eta[i][j], beta);
    double sum = 0.0;
    for (int k = 0; k < V; k++)
        if (!tabu[k])
            sum += pow(tau[i][k], alpha) * pow(eta[i][k], beta);
    P /= sum;
    return P;
}

void upd_pheromone(struct Path *P, int V, int Q, double p, double ***tau){
    double dlt;
    for (int i = 0; i < V; i++)
        for (int j = 0; j < V; j++)
            if (i != j)
                (*tau)[i][j] = (*tau)[i][j] * (1.0 - p);
    for (int k = 0; k < V; k++)
        for (int i = 0; i < V; i++){
            dlt = (double) Q / (double) P[k].l;
            (*tau)[P[k].A[i]][P[k].A[(i + 1) % V]] += dlt;
            (*tau)[P[k].A[(i + 1) % V]][P[k].A[i]] += dlt;
        }
}

struct Path ACO_solve(struct Graph *G, int V, int alpha, int beta, double p, double tau0, double rf, int t){
    int Q = 0, i, j, *tabu, rv;
    double sum_p, r;
    struct Path BP;
    BP.A = malloc(V * sizeof(int));
    for (i = 0; i < V; i++){
        BP.A[i] = i;
        Q += G -> A[i][(i + 1) % V];
    }
    BP.l = Q;
    struct Path *P = malloc(V * sizeof(struct Path));
    for (int k = 0; k < V; k++){
        P[k].A = malloc(V * sizeof(int));
        P[k].l = 0;
    }
    double **tau = (double **) malloc(V * sizeof(double *));
    double **eta = (double **) malloc(V * sizeof(double *));
    for (i = 0; i < V; i++){
        tau[i] = calloc(V, sizeof(double));
        eta[i] = calloc(V, sizeof(double));
    }
    for (i = 0; i < V; i++)
        for (j = 0; j < V; j++){
            if (i != j){
                tau[i][j] = tau0;
                eta[i][j] = 1 / (double) (G -> A[i][j]);
            }
        }
    for (t; t; t--){
        for (int k = 0; k < V; k++){
            i = 0;
            tabu = (int *) calloc(V, sizeof(int));
            P[k].A[i] = k;
            tabu[k]++;
            while (i < V){
                sum_p = 0.0;
                rv = rand() % V;
                if (((double) rand() / (double) RAND_MAX < rf) && !tabu[rv]){
                    P[k].l += G -> A[P -> A[i]][rv];
                    P[k].A[++i] = rv;
                    tabu[rv]++;
                }
                else {
                    j = 0;
                    r = (double) rand() / (double) RAND_MAX;
                    for (j; j < V && sum_p < r; j++)
                        if (!tabu[j])
                            sum_p += probability(P[k].A[i], j, tabu, V, alpha, beta, tau, eta);
                    if (sum_p)
                        j = (j - 1) % V;
                    else
                        j = k;
                    P[k].l += G -> A[P[k].A[i]][j];
                    P[k].A[++i] = j;
                    tabu[j]++;
                }
            }
            if (P[k].l < BP.l){
                BP.l = P[k].l;
                for (i = 0; i < V; i++)
                    BP.A[i] = P[k].A[i];
            }
            free(tabu);
        }
        upd_pheromone(P, V, Q, p, &tau);
    }
    for (i = 0; i < V; i++){
        free(tau[i]);
        free(eta[i]);
    }
    free(P);
    free(tau);
    free(eta);
    return BP;
}

int main(void){
    int V, E, src, dst, wt, alpha, beta, t;
    double p, tau0, rf;
    printf("Number of vertecies: ");
    scanf("%i", &V);
    struct Graph *G = create_graph(V);
    printf("Number of edges: ");
    scanf("%i", &E);
    for (E; E; E--){
        printf("Enter the src vertex, dst vertex and edge weight: ");
        scanf("%i %i %i", &src, &dst, &wt);
        create_edge(G, src, dst, wt);
    }
    printf("\nEntered graph:\n");
    print_graph(G, V);
    printf("\nEnter ACO algorithm parameters:\n");
    printf("alpha (pheromone wt) and beta (visibility wt): ");
    scanf("%i %i", &alpha, &beta);
    printf("p (pheromone valitilization coeff.) and tau0 (initial pheromone track): ");
    scanf("%lf %lf", &p, &tau0);
    printf("rf (random factor for movement): ");
    scanf("%lf", &rf);
    printf("Number of iterations: ");
    scanf("%i", &t);
    struct Path P = ACO_solve(G, V, alpha, beta, p, tau0, rf, t);
    output_path(&P, V);
    delete_graph(G, V);
    system("Pause");
    return 0;
}