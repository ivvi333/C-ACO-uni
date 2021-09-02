#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <time.h>
#include "../Graphs/graph.h"

// Стуктура путь
struct Path{
    int l; // общая длина пути
    int *A; // массив вершин пути
};

void delete_path_arr(struct Path *P, int V){
    for (register int i = 0; i < V; i++)
        free(P[i].A);
    free(P);
}

void output_path(struct Path *P, int V){
    printf("\nPath: %i ", P -> A[0]);
    for (register int i = 1; i <= V; i++)
        printf("--> %i ", P -> A[i]);
    printf("\nPath length: %i\n", P -> l);
}

double probability(int i, int j, int *tabu, int V, int alpha, int beta, double **tau, double **eta){
    double P = pow(tau[i][j], (double) alpha) * pow(eta[i][j], (double) beta);
    double sum = 0.0;
    for (register int k = 0; k < V; k++)
        if (!tabu[k])
            sum += pow(tau[i][k], alpha) * pow(eta[i][k], beta);
    P /= sum;
    return P;
}

void upd_pheromone(struct Path *P, int V, int Q, double p, double ***tau){
    double dlt;
    for (register int i = 0; i < V; i++)
        for (register int j = 0; j < V; j++)
            (*tau)[i][j] = (*tau)[i][j] * (1.0 - p);
    for (register int k = 0; k < V; k++){
        dlt = (double) Q / (double) P[k].l;
        for (register int i = 0; i < V; i++){
            (*tau)[P[k].A[i]][P[k].A[i + 1]] += dlt;
            (*tau)[P[k].A[i + 1]][P[k].A[i]] += dlt;
        }
    }
}

struct Path ACO_solve(struct Graph *G, int V, int alpha, int beta, double p, double tau0, int t){
    int Q = 90;
    register int i, j, k;
    double sum_p, r;
    struct Path BP;
    BP.A = malloc((V + 1) * sizeof(int));
    BP.l = 0;
    for (i = 0; i < V; i++){
        BP.A[i] = i;
        BP.l += G -> A[i][(i + 1) % V];
    }
    BP.A[i] = 0;
    struct Path *P = malloc(V * sizeof(struct Path));
    for (k = 0; k < V; k++){
        P[k].A = malloc((V + 1) * sizeof(int));
        P[k].l = 0;
    }
    int *tabu = (int *) calloc(V, sizeof(int));
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
                eta[i][j] = 1.0 / (double) (G -> A[i][j]);
            }
        }
    for (t; t; t--){
        for (k = 0; k < V; k++){
            i = 0;
            P[k].l = 0;
            P[k].A[i++] = k;
            tabu[k]++;
            while (i < V){
                sum_p = 0.0;
                r = (double) rand() / (double) RAND_MAX;
                while (r == 0.0 || r == 1.0)
                    r = (double) rand() / (double) RAND_MAX;
                for (j = 0; j < V && sum_p < r; j++)
                    if (!tabu[j])
                        sum_p += probability(P[k].A[i - 1], j, tabu, V, alpha, beta, tau, eta);
                j = (j - 1) % V;
                P[k].l += G -> A[P[k].A[i - 1]][j];
                P[k].A[i++] = j;
                tabu[j]++;
            }
            P[k].l += G -> A[P[k].A[i - 1]][k];
            P[k].A[i] = k;
            if (P[k].l < BP.l){
                BP.l = P[k].l;
                memmove(BP.A, P[k].A, (V + 1) * sizeof(int));
            }
            memset(tabu, 0, V * sizeof(int));
        }
        upd_pheromone(P, V, Q, p, &tau);
    }
    for (i = 0; i < V; i++){
        free(tau[i]);
        free(eta[i]);
    }
    delete_path_arr(P, V);
    free(tau);
    free(eta);
    free(tabu);
    return BP;
}

int main(void){
    srand(time(NULL));
    int V, alpha, beta, t;
    double p, tau0;
    printf("Number of vertecies: ");
    scanf("%i", &V);
    struct Graph *G = create_graph(V);
    input_graph(G);
    printf("\nEnter ACO algorithm parameters:\n");
    printf("alpha (pheromone wt) and beta (visibility wt): ");
    scanf("%i %i", &alpha, &beta);
    printf("p (pheromone valitilization coeff.) and tau0 (initial pheromone track): ");
    scanf("%lf %lf", &p, &tau0);
    printf("Number of iterations: ");
    scanf("%i", &t);
    clock_t begin = clock();
    struct Path P = ACO_solve(G, V, alpha, beta, p, tau0, t);
    clock_t end = clock();
    output_path(&P, V);
    printf("Execution time: %lf s\n", (double)(end - begin) / CLOCKS_PER_SEC);
    free(P.A);
    delete_graph(G);
    system("Pause");
    return 0;
}