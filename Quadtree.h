/*
Copyright (c) 2013, Deepak Narayanan & Quan Ngyugen
*/

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

#define N 60
#define INSERT_COARSE_LIM 120

// Definition of a line node, which is fundamental unit of our linked list
struct line_node {
  struct line_node* next;
  Line* line;
};
typedef struct line_node line_node;

line_node* line_node_new(Line* line);

struct line_list {
  size_t num_lines;
  struct line_node* head;
  struct line_node* tail;
};
typedef struct line_list line_list;

struct quad_tree {
  struct quad_tree *quad1, *quad2, *quad3, *quad4;
  line_node* lines;
  size_t num_lines;  // total lines contained, not the length of 'lines'.
  double xmin, xmax, ymin, ymax;
};
typedef struct quad_tree quad_tree;

line_list *line_list_new();
quad_tree *quad_tree_new(double xmin, double xmax, double ymin, double ymax);

// void line_list_delete(line_list* list);
// void quad_tree_delete(quad_tree * tree);

// Inserts a new line into the given linked list, making sure that
// the input line is not modified by this operation in any way
void insert_line(line_node** lines, line_node* new_line);
// Merges the two lists, does not modify list2
void merge_lists(line_node** list1, line_node* list2);

int get_quad_type_line(Vec p1, Vec p2, quad_tree* tree);
int get_quad_type(quad_tree* tree, line_node* node, double timeStep);

void quadtree_insert_lines(quad_tree* tree, line_node* new_lines, double timeStep, int num_lines);

#endif
