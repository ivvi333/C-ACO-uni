#include <stdio.h>
#include <stdlib.h>
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

struct Path *ACO_solve(struct Graph *G, int V, int alpha, int beta, double p, double tau0, int t){
    int Q = 0;
    struct Path *P = malloc(sizeof(struct Path));
    P -> V = V;
    P -> A = malloc(V * sizeof(int));
    for (int i = 0; i < V; i++){
        Q += G -> A[i][(i + 1) % V];
    }
    P -> l = Q;
    double **tau = (double **) malloc(V * sizeof(double *));
    double **eta = (double **) malloc(V * sizeof(double *));
    for (int i = 0; i < V; i++){
        tau[i] = calloc(V, sizeof(double));
        eta[i] = calloc(V, sizeof(double));
    }
    for (int i = 0; i < V; i++)
        for (int j = 0; j < V; j++){
            if (i != j){
                tau[i][j] = tau0;
                eta[i][j] = 1/(G -> A[i][j]);
            }
        }
    
    for (int i = 0; i < V; i++){
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