/*
 * user_types.hpp
 *
 *  Created on: 11 Oct 2021
 *      Author: mndx
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

#include <vector>

const int inf = 3e+8;
extern int tot_num_ops;

typedef struct FibHeapProperties {
    bool deg_is_num_child;
    int num_nodes;
} fib_props;

typedef struct Node {
    Node* left;
    Node* right;
    Node* p;
    Node* child;

    std::vector<int> adj_nodes;

    int key;
    int degree;
    int index;
    bool mark;
} node;

class FibHeap {
public:
    int n;
    node* min;
    FibHeap() { min = NULL; n = 0; }
};

#endif /* USER_TYPES_HPP_ */
