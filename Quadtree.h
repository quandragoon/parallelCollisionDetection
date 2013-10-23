/**
 * Quadtree.c -- Quadtree data structure which contains line segments
 * Copyright (c) 2012 the Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 **/

#ifndef QUADTREE_H_
#define QUADTREE_H_

#include <stdlib.h>
#include "./Line.h"
#include "./Vec.h"

#define MUL_TYPE 0
#define Q1_TYPE 1
#define Q2_TYPE 2
#define Q3_TYPE 3
#define Q4_TYPE 4

#define N 62
#define INSERT_COARSE_LIM 120

// Definition of a line node, which is fundamental unit of our linked list
struct line_node {
  struct line_node* next;
  Line* line;
};
typedef struct line_node line_node;

line_node* line_node_new(Line* line);

struct quad_tree {
  // quad1 is top left, quad2 is top right, quad3 is bottom left and
  // quad4 is bottom right.
  struct quad_tree *quad1, *quad2, *quad3, *quad4;
  // Head of linked list that contains all line segments stored at that level
  // of the tree
  line_node* lines;
  size_t num_lines;  // total lines contained, not the length of 'lines'.
  // Coordinates of bounding box of the quadtree
  double xmin, xmax, ymin, ymax;
};
typedef struct quad_tree quad_tree;

quad_tree *quad_tree_new(double xmin, double xmax, double ymin, double ymax);

void quad_tree_delete(quad_tree * tree);

// Inserts a new line into the given linked list, making sure that
// the input line is not modified by this operation in any way
void insert_line(line_node** lines, line_node* new_line);
// Merges the two lists, does not modify list2
void merge_lists(line_node** list1, line_node* list2);

int get_quad_type_line(Vec p1, Vec p2, quad_tree* tree);
int get_quad_type(quad_tree* tree, line_node* node, double timeStep);

void quadtree_insert_lines(quad_tree* tree, line_node* new_lines, double timeStep, int num_lines);

#endif  // QUADTREE_H_
