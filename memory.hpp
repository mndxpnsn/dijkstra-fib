/*
 * memory.hpp
 *
 *  Created on: 11 Oct 2021
 *      Author: mndx
 */

#ifndef MEMORY_HPP_
#define MEMORY_HPP_

#include "user_types.hpp"

bool** bool2D(const int size);
int** int2D(const int size);
void free_bool2D(bool** p, int size) ;
void free_int2D(int** p, int size);
void free_node_ref(node** v_ref, int size);

#endif /* MEMORY_HPP_ */
