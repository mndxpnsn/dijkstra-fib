/*
 * memory.cpp
 *
 *  Created on: 11 Oct 2021
 *      Author: mndx
 */

#include "user_types.hpp"

bool** bool2D(const int size) {
    bool** p = new bool*[size];

    for(int i = 0; i < size; ++i) {
        p[i] = new bool[size];
    }

    return p;
}

int** int2D(const int size) {
    int** p = new int*[size];

    for(int i = 0; i < size; ++i) {
        p[i] = new int[size];
    }

    return p;
}

void free_bool2D(bool** p, int size) {
    for(int i = 0; i < size; ++i) {
        delete [] p[i];
    }

    delete [] p;
}

void free_int2D(int** p, int size) {
    for(int i = 0; i < size; ++i) {
        delete [] p[i];
    }

    delete [] p;
}

void free_node_ref(node** v_ref, int size) {
    for(int i = 0; i < size; ++i) {
        delete v_ref[i];
    }

    delete [] v_ref;
}



