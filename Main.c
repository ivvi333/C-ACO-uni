#include <stdio.h>
#include "../Graphs/graph.h"

int main(void){
    int V;
    printf("Number of vertecies: ");
    scanf("%i", &V);
    struct Graph *G = create_graph(V);
    return 0;
}