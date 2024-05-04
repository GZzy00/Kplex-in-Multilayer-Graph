#include "Graph.h"
#include "Utility.h"
#include "Timer.h"

using namespace std;

int main(int argc, char *argv[])
{
    int K = 5, Sigma = 2;
    Graph *graph = new Graph("./datasets/homo.txt", K, Sigma);
    graph->read_graph();
    graph->kPlex_exact();
    // for (int k = 2; k < 6; k ++)
    // {
    //     for (int s = 2; s < 4; s++)
    //     {
    //         K = k, Sigma = s;
    //         Graph *graph = new Graph("./datasets/slashdot.txt", K, Sigma);
    //         graph->read_graph();
    //         graph->kPlex_exact();
    //     }
    // }
    return 0;
}