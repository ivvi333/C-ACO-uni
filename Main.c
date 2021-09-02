#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <time.h>
#include "../Graphs/graph.h"

// Стуктура муравей
struct Ant{
    int l; // общая длина пути муравья
    int *P; // массив вершин пути муравья
};

void delete_ant_arr(struct Ant *A, int V){
    for (register int i = 0; i < V; i++)
        free(A[i].P);
    free(A);
}

void output_path(struct Ant *A, int V){
    printf("\nPath: %i ", A -> P[0]);
    for (register int i = 1; i <= V; i++)
        printf("--> %i ", A -> P[i]);
    printf("\nPath length: %i\n", A -> l);
}

double probability(int i, int j, int *tabu, int V, int alpha, int beta, double **tau, double **eta){
    double pr = pow(tau[i][j], (double) alpha) * pow(eta[i][j], (double) beta);
    double sum = 0.0;
    for (register int k = 0; k < V; k++)
        if (!tabu[k])
            sum += pow(tau[i][k], alpha) * pow(eta[i][k], beta);
    pr /= sum;
    return pr;
}

void upd_pheromone(struct Ant *A, int V, int Q, double p, double ***tau){
    double dlt;
    for (register int i = 0; i < V; i++)
        for (register int j = 0; j < V; j++)
            (*tau)[i][j] = (*tau)[i][j] * (1.0 - p);
    for (register int k = 0; k < V; k++){
        dlt = (double) Q / (double) A[k].l;
        for (register int i = 0; i < V; i++){
            (*tau)[A[k].P[i]][A[k].P[i + 1]] += dlt;
            (*tau)[A[k].P[i + 1]][A[k].P[i]] += dlt;
        }
    }
}

struct Ant ACO_solve(struct Graph *G, int V, int alpha, int beta, double p, double tau0, int t){
    const int Q = 90;
    register int i, j, k;
    double sum_p, r;
    struct Ant BA;
    BA.P = malloc((V + 1) * sizeof(int));
    BA.l = 0;
    for (i = 0; i < V; i++){
        BA.P[i] = i;
        BA.l += G -> A[i][(i + 1) % V];
    }
    BA.P[i] = 0;
    struct Ant *A = malloc(V * sizeof(struct Ant));
    for (k = 0; k < V; k++){
        A[k].P = malloc((V + 1) * sizeof(int));
        A[k].l = 0;
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
            A[k].l = 0;
            A[k].P[i++] = k;
            tabu[k]++;
            while (i < V){
                sum_p = 0.0;
                r = (double) rand() / (double) RAND_MAX;
                while (r == 0.0 || r == 1.0)
                    r = (double) rand() / (double) RAND_MAX;
                for (j = 0; j < V && sum_p < r; j++)
                    if (!tabu[j])
                        sum_p += probability(A[k].P[i - 1], j, tabu, V, alpha, beta, tau, eta);
                j = (j - 1) % V;
                A[k].l += G -> A[A[k].P[i - 1]][j];
                A[k].P[i++] = j;
                tabu[j]++;
            }
            A[k].l += G -> A[A[k].P[i - 1]][k];
            A[k].P[i] = k;
            if (A[k].l < BA.l){
                BA.l = A[k].l;
                memmove(BA.P, A[k].P, (V + 1) * sizeof(int));
            }
            memset(tabu, 0, V * sizeof(int));
        }
        upd_pheromone(A, V, Q, p, &tau);
    }
    for (i = 0; i < V; i++){
        free(tau[i]);
        free(eta[i]);
    }
    delete_ant_arr(A, V);
    free(tau);
    free(eta);
    free(tabu);
    return BA;
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
    struct Ant A = ACO_solve(G, V, alpha, beta, p, tau0, t);
    clock_t end = clock();
    output_path(&A, V);
    printf("Execution time: %lf s\n", (double)(end - begin) / CLOCKS_PER_SEC);
    free(A.P);
    delete_graph(G);
    system("Pause");
    return 0;
}