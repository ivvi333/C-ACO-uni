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

// Удаляет массив структур муравей по адресу A размера V
void delete_ant_arr(struct Ant *A, int V){
    for (register int i = 0; i < V; i++)
        free(A[i].P);
    free(A);
}

// Выводит маршрут муравья по адресу A;
// V + 1 - кол-во посещённых муравьем вершин 
void output_path(struct Ant *A, int V){
    printf("\nPath: %i ", A -> P[0]);
    for (register int i = 1; i <= V; i++)
        printf("--> %i ", A -> P[i]);
    printf("\nPath length: %i\n", A -> l);
}

// Рассчитывает вероятность перехода муравья из города i в город j
double probability(int i, int j, int *tabu, int V, int alpha, int beta, double **tau, double **eta){
    double pr = pow(tau[i][j], (double) alpha) * pow(eta[i][j], (double) beta);
    double sum = 0.0;
    for (register int k = 0; k < V; k++)
        if (!tabu[k])
            sum += pow(tau[i][k], alpha) * pow(eta[i][k], beta);
    pr /= sum;
    return pr;
}

// Обновляет виртуальный феромонный след
void upd_pheromone(struct Ant BA, int V, double rho, double tau_min, double tau_max, double **tau){
    double dlt;
    for (register int i = 0; i < V; i++)
        for (register int j = 0; j < V; j++){
            tau[i][j] = tau[i][j] * rho;
            if (tau[i][j] < tau_min)
                tau[i][j] = tau_min;
        }
    dlt = 1.0 / (double) BA.l;
    for (register int i = 0; i < V; i++){
        tau[BA.P[i]][BA.P[i + 1]] += dlt;
        tau[BA.P[i + 1]][BA.P[i]] += dlt;
        if (tau[BA.P[i]][BA.P[i + 1]] > tau_max){
            tau[BA.P[i]][BA.P[i + 1]] = tau_max;
            tau[BA.P[i + 1]][BA.P[i]] = tau_max;
        }
    }
}

struct Ant ACO_solve(struct Graph *G, int V, int alpha, int beta, double rho, double p_best, int t){
    // Выделение памяти, инициализация и определение переменных
    const int T = t;
    register int i, j, k;
    double tau_max, tau_min, rt_p_best, sum_p, r;
    // Создаём глобально лучшего и итерационно лучшего муравья 
    struct Ant BG, BI;
    BG.P = malloc((V + 1) * sizeof(int));
    BI.P = malloc((V + 1) * sizeof(int));
    BG.l = 0;
    // В начале алгоритма глобально лучшим будем считать маршрут 0 -> 1 -> 2 -> ... -> V - 1 -> 0
    for (i = 0; i < V; i++){
        BG.P[i] = i;
        BG.l += G -> A[i][(i + 1) % V];
    }
    BG.P[i] = 0;
    // Определяем начальные границы
    tau_max = 1.0 / (1.0 - rho) * 1.0 / (double) BG.l;
    rt_p_best = pow(p_best, 1.0 / (double) V);
    tau_min = tau_max * (1 - rt_p_best) / (((double) V / 2.0) * rt_p_best);
    // A - массив муравьёв
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
                tau[i][j] = tau_max;
                eta[i][j] = 1.0 / (double) (G -> A[i][j]);
            }
        }

    // Основной цикл, выполняется t раз
    for (t; t; t--){
        // В начале каждой итерации длина итерационно лучшего маршрута максимальна, т.е. наихудшая
        BI.l = __INT_MAX__;
        // Для каждого муравья в городе k
        for (k = 0; k < V; k++){
            i = 0;
            A[k].l = 0;
            A[k].P[i++] = k;
            tabu[k]++;
            // Пока маршрут, за исключением последней вершины, не получен
            while (i < V){
                sum_p = 0.0;
                r = (double) rand() / (double) RAND_MAX;
                while (r == 0.0 || r == 1.0)
                    r = (double) rand() / (double) RAND_MAX;
                // Пока есть города или пока сумма вероятностей меньше случайного числа r - "принцип рулетки"
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
            // Если полученный маршрут короче итерационно лучшего
            if (A[k].l < BI.l){
                BI.l = A[k].l;
                memmove(BI.P, A[k].P, (V + 1) * sizeof(int));
                // Если итерационно лучший маршрут короче глобально лучшего
                if (BI.l < BG.l){
                    BG.l = BI.l;
                    memmove(BG.P, BI.P, (V + 1) * sizeof(int));
                    // Обновляем границы
                    tau_max = 1.0 / (1.0 - rho) * 1.0 / (double) BG.l;
                    rt_p_best = pow(p_best, 1.0 / (double) V);
                    tau_min = tau_max * (1 - rt_p_best) / (((double) V / 2.0) * rt_p_best);
                }
            }
            // Очищаем список табу
            memset(tabu, 0, V * sizeof(int));
        }
        // Обновляем феромон
        upd_pheromone(BI, V, rho, tau_min, tau_max, tau);
    }

    // Освобождение памяти
    for (i = 0; i < V; i++){
        free(tau[i]);
        free(eta[i]);
    }
    delete_ant_arr(A, V);
    free(tau);
    free(eta);
    free(tabu);
    free(BI.P);
    return BG;
}

int main(void){
    srand(time(NULL));
    FILE *fp = NULL;
    FILE *gnuplot = NULL;
    int V, alpha, beta, t, tests;
    double rho, p_best;
    printf("Number of vertecies: ");
    scanf("%i", &V);
    struct Graph *G = create_graph(V);
    input_graph(G);
    printf("\nEnter ACO algorithm parameters:\n");
    printf("alpha (pheromone wt) and beta (visibility wt): ");
    scanf("%i %i", &alpha, &beta);
    printf("rho (pheromone valitilization coeff.) and p_best (expected best solution probability): ");
    scanf("%lf %lf", &rho, &p_best);
    printf("Number of iterations: ");
    scanf("%i", &t);
    printf("Number of tests: ");
    scanf("%i", &tests);
    fp = fopen("data.tmp", "w");
    gnuplot = _popen("gnuplot -persistent", "w");
    char *GnuCommands[] = {"set title \"Results\"", "plot 'data.tmp'"};
    for (int test_num = 0; test_num < tests; test_num++){
        clock_t begin = clock();
        struct Ant A = ACO_solve(G, V, alpha, beta, rho, p_best, t);
        fprintf(fp, "%d %d\n", test_num, A.l);
        clock_t end = clock();
        output_path(&A, V);
        printf("Execution time: %lf s\n", (double)(end - begin) / CLOCKS_PER_SEC);
        free(A.P);
    }
    for (register int i = 0; i < 2; i++)
        fprintf(gnuplot, "%s\n", GnuCommands[i]);
    delete_graph(G);
    fclose(fp);
    _pclose(gnuplot);
    system("Pause");
    return 0;
}