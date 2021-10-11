/*
 * fib_heap_support.hpp
 *
 *  Created on: 11 Oct 2021
 *      Author: mndx
 */

#ifndef FIB_HEAP_SUPPORT_HPP_
#define FIB_HEAP_SUPPORT_HPP_

#include "user_types.hpp"

void print_root_list(node* z);
void print_child_list(node* child);
void print_list(node* z);
bool check_fib_heap(FibHeap* H);

#endif /* FIB_HEAP_SUPPORT_HPP_ */
