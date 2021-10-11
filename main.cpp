/*
 * main.cpp
 *
 *  Created on: 11 Oct 2021
 *      Author: mndx
 */

#include <iostream>
#include <math.h>
#include <vector>

#include "fib_heap.hpp"
#include "user_types.hpp"

int main(int argc, char* argv[]) {

    //Declarations
    int s = 2; //Start vertex. The minimum vertex has index 1
    int n = 2499; //Number of vertices
    int num_edges = 3125; //Number of edges

    //Create edges
    std::vector< std::vector<int> > edges;
    for(int i = 0; i < num_edges; ++i) {
        int start_vert = rand() % n + 1;
        int end_vert = rand() % n + 1;
        int weight = rand() % 200 + 1;

        std::vector<int> edge_elem;
        edge_elem.push_back(start_vert);
        edge_elem.push_back(end_vert);
        edge_elem.push_back(weight);
        edges.push_back(edge_elem);
    }

    //Compute distances to nodes from start vertex
    std::vector<int> results = shortest_reach(n, edges, s);

    //Print results
    float tot_num_ops_est = 5*n + 4*num_edges + 6.4*n*log(n)/log(2);
    float complexity_ratio = tot_num_ops / tot_num_ops_est;
    int size_results = (int) results.size();
    for(int i = 0; i < size_results; ++i) {
        std::cout << results[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "number of operations estimated 5V + 4E + 6.4VlgV: " << tot_num_ops_est << std::endl;
    std::cout << "number of operations measured: " << tot_num_ops << std::endl;
    std::cout << "complexity ratio: " << complexity_ratio << std::endl;
    std::cout << "done" << std::endl;

    return 0;
}
