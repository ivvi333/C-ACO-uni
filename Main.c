#include <stdio.h>
#include <windows.h>
#include "../Graphs/graph.h"

int main(void){
    int V, E, src, dst, wt, alpha, beta, Q, tau0;
    double p;
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
    printf("Q (optimal path length): ");
    scanf("%i", &Q);
    printf("p (pheromone valitilization coeff.) and tau0 (initial pheromone track): ");
    scanf("%lf %i", &p, &tau0);
    delete_graph(G, V);
    system("Pause");
    return 0;
}