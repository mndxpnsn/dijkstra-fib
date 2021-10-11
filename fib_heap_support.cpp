/*
 * fib_heap_support.cpp
 *
 *  Created on: 11 Oct 2021
 *      Author: mndx
 */

#include <iostream>

#include "user_types.hpp"

void print_root_list(node* z) {
    node* xt = z;
    if(xt != NULL) {
        do {
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
            xt = xt->right;
        } while(xt != z);
    }
}

void print_child_list(node* child) {
    node* xt = child;
    if(xt != NULL) {
        do {
            std::cout << "xt->child->key: " << xt->key;
            std::cout << ", xt->child->degree: " << xt->degree << std::endl;
            if(xt->child != NULL) {
                std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
            }
            xt = xt->right;
        } while(xt != child);
    }
}

void print_list(node* z) {
    node* xt = z;
    if(xt != NULL) {
        do {
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
            if(xt->child != NULL) {
                print_child_list(xt->child);
            }
            xt = xt->right;
        } while(xt != z);
    }
}

bool numbers_children_match(node* z, int& num_nodes) {
    bool nums_match = true;
    int num_of_nodes = 0;

    node* xt = z->child;
    if(xt != NULL) {
        do {
            num_of_nodes++;
            if(xt->child != NULL) {
                nums_match = numbers_children_match(xt, num_nodes);
                if(!nums_match) { return false; }
            }
            xt = xt->right;
        } while(xt != z->child);

        num_nodes = num_nodes + num_of_nodes;

        if(num_of_nodes == z->degree) { nums_match = true; }
        else { nums_match = false; }
    }

    return nums_match;
}

fib_props numbers_match(node* z) {

    bool nums_match = true;
    int num_nodes = 0;
    fib_props fib_heap_props = { nums_match, num_nodes };

    node* xt = z;
    if(xt != NULL) {
        do {
            num_nodes++;
            nums_match = numbers_children_match(xt, num_nodes);
            fib_heap_props.deg_is_num_child = nums_match;
            fib_heap_props.num_nodes = num_nodes;
            if(!nums_match) { return fib_heap_props; }
            xt = xt->right;
        } while(xt != z);
    }

    fib_heap_props.deg_is_num_child = nums_match;
    fib_heap_props.num_nodes = num_nodes;

    return fib_heap_props;
}

bool is_fib_heap_children(node* z) {
    bool is_fibheap = true;

    node* xt = z->child;
    if(xt != NULL) {
        do {
            if(xt->p->key > xt->key) {
                return is_fibheap = false;
            }
            if(xt->child != NULL) {
                is_fibheap = is_fib_heap_children(xt);
                if(!is_fibheap) { return false; }
            }
            xt = xt->right;
        } while(xt != z->child);
    }

    return is_fibheap;
}

bool is_fib_heap(node* z) {
    bool is_fibheap = true;

    node* xt = z;
    if(xt != NULL) {
        do {
            is_fibheap = is_fib_heap_children(xt);
            if(!is_fibheap) { return false; }
            xt = xt->right;
        } while(xt != z);
    }

    return is_fibheap;
}

bool check_fib_heap(FibHeap* H) {
    //This is the general test for the fibonacci heap.
    //The function returns true if the heap satisfies
    //the fibonacci heap properties

    //Compute heap properties
    fib_props fh_props = numbers_match(H->min);
    bool heap_is_fibheap = is_fib_heap(H->min);

    //Check if number of children equal node degrees
    bool deg_is_num_childs = fh_props.deg_is_num_child;

    //Check if number of nodes counted in heap equals H.n
    int num_nodes = fh_props.num_nodes;
    bool num_nodes_match = (num_nodes == H->n);

    //Check to see if heap is properly structured
    bool heap_is_ok = num_nodes_match && deg_is_num_childs && heap_is_fibheap;

    return heap_is_ok;
}
