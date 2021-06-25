#include <stdio.h>
#include "../Graphs/graph.h"

int main(void){
    int V, E, src, dst, wt;
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
    print_graph(G, V);
    delete_graph(G, V);
    return 0;
}