#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include "../Graphs/graph.h"

// Стуктура путь
struct Path{
    int V; // количество вершин
    int l; // общая длина пути
    int *A; // массив вершин пути
};

void delete_path(struct Path *P){
    free(P -> A);
    free(P);
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

struct Path *ACO_solve(struct Graph *G, int V, int alpha, int beta, double p, double tau0, int t){
    int Q = 0, is_atStart, i, j, *tabu;
    struct Path *P = malloc(sizeof(struct Path));
    P -> V = V;
    P -> A = malloc(V * sizeof(int));
    for (i = 0; i < V; i++){
        Q += G -> A[i][(i + 1) % V];
    }
    P -> l = Q;
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
                eta[i][j] = 1/(G -> A[i][j]);
            }
        }
    for (t; t; t--){
        for (int k = 0; k < V; k++){
            i = 0;
            tabu = (int *) calloc(V, sizeof(int));
            tabu[k]++;
            is_atStart = 0;
            while (!is_atStart){
                j = 0;
                for (j; j < V; j++){
                    if (!tabu[j]){
                        
                    }
                }
            }
        }
    }
    for (i = 0; i < V; i++){
        free(tau[i]);
        free(eta[i]);
    }
    free(tau);
    free(eta);
    return P;
}

int main(void){
    int V, E, src, dst, wt, alpha, beta, t;
    double p, tau0;
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
    printf("Number of iterations: ");
    scanf("%i", &t);
    struct Path *P = ACO_solve(G, V, alpha, beta, p, tau0, t);
    delete_graph(G, V);
    delete_path(P);
    system("Pause");
    return 0;
}