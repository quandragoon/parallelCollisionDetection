#include "./Quadtree.h"
#include "./Line.h"
#include "./Vec.h"
#include <stdlib.h>
#include <stdio.h>

#include "assert.h"

line_node* line_node_new(Line* line) {
  line_node *const new_line = (line_node *const)malloc(sizeof(struct line_node));
  new_line->line = line;
  return new_line;
}

quad_tree *quad_tree_new(double xmin, double xmax, double ymin, double ymax) {
  quad_tree * root = (quad_tree *)malloc(sizeof(quad_tree));
  root->quad1 = root->quad2 = root->quad3 = root->quad4 = NULL;
  root->lines = NULL;
  root->num_lines = 0;
  root->xmin = xmin;
  root->xmax = xmax;
  root->ymin = ymin;
  root->ymax = ymax;
  return root;
}

/*void quad_tree_delete(quad_tree * tree) {
  //Call recursively
  if (tree->quad1 != NULL)
    quad_tree_delete(tree->quad1);
  if (tree->quad2 != NULL)
    quad_tree_delete(tree->quad2);
  if (tree->quad3 != NULL)
    quad_tree_delete(tree->quad3);
  if (tree->quad4 != NULL)
    quad_tree_delete(tree->quad4);
  free(tree->lines);
}*/

// Inserts a new line into the given linked list, making sure that
// the input line is not modified by this operation in any way
void insert_line(line_node** lines, line_node* new_line) {
  if (*lines == NULL) {
    new_line->next = NULL;
    *lines = new_line;
  } else {
    new_line->next = *lines;
    *lines = new_line;
  }
}

// Merges the two lists, does not modify list2
void merge_lists(line_node** list1, line_node* list2) {
  if (list2 == NULL) return;
  if (*list1 == NULL) {
    *list1 = list2;
  } else {
    line_node* node = *list1;
    line_node* prev = NULL;
    while (node != NULL) {
      prev = node;
      node = node->next;
    }
    if (prev != NULL)
      prev->next = list2;
  }
}

int get_quad_type_line(Vec p1, Vec p2, quad_tree* tree) {
  double xmin = tree->xmin;
  double xmax = tree->xmax;
  double ymin = tree->ymin;
  double ymax = tree->ymax;
  double xmid = (xmin + xmax) / 2.0;
  double ymid = (ymin + ymax) / 2.0;

  if (!(((p1.x - xmid)*(p2.x - xmid) > 0) && ((p1.y - ymid)*(p2.y - ymid) > 0)))
    return MUL_TYPE;

  int xid = (p1.x - xmid > 0) ? 1 : 0;
  int yid = (p1.y - ymid > 0) ? 1 : 0;
  int quad = 2 * yid + xid + 1;
  return quad;
}

int get_quad_type(quad_tree* tree, line_node* node, double timeStep) {
  Line* actual_line = node->line;
  Vec p1 = actual_line->p1;
  Vec p2 = actual_line->p2;

  Vec new_p1 = Vec_add(p1, Vec_multiply(node->line->velocity, timeStep));
  Vec new_p2 = Vec_add(p2, Vec_multiply(node->line->velocity, timeStep));

  int first_quad = get_quad_type_line(p1, p2, tree);
  int second_quad = get_quad_type_line(new_p1, new_p2, tree);
  return (first_quad == second_quad) ? first_quad : MUL_TYPE;
}

void quadtree_insert_lines(quad_tree* tree, line_node* new_lines, double timeStep, int num_lines) {
  tree->num_lines = num_lines;
  double xmax = tree->xmax;
  double xmin = tree->xmin;
  double ymax = tree->ymax;
  double ymin = tree->ymin;

  if (num_lines <= N) {
    tree->lines = new_lines;
    return;
  }

  line_node *quad1, *quad2, *quad3, *quad4, *lines;
  quad1 = quad2 = quad3 = quad4 = lines = NULL;
  int num_quad1, num_quad2, num_quad3, num_quad4, num_parent_lines;
  num_quad1 = num_quad2 = num_quad3 = num_quad4 = num_lines = 0;

  line_node* cur = new_lines;
  line_node* next;
  int type;
  while (cur != NULL) {
    type = get_quad_type(tree, cur, timeStep);
    next = cur->next;
    switch (type) {
      case Q1_TYPE:
        insert_line(&quad1, cur);
        num_quad1++;
        break;
      case Q2_TYPE:
        insert_line(&quad2, cur);
        num_quad2++;
        break;
      case Q3_TYPE:
        insert_line(&quad3, cur);
        num_quad3++;
        break;
      case Q4_TYPE:
        insert_line(&quad4, cur);
        num_quad4++;
        break;
      case MUL_TYPE:
        insert_line(&lines, cur);
        num_parent_lines++;
        break;
      default:
        return;
    }
    cur = next;
  }

  double xmid = (xmin + xmax) / 2.0;
  double ymid = (ymin + ymax) / 2.0;
  tree->lines = lines;
  /*
  if (quad1->head) {
    tree->quad1 = quad_tree_new(xmin, xmid, ymin, ymid);
  if (quad1->num_lines > INSERT_COARSE_LIM)
    cilk_spawn quadtree_insert_lines(tree->quad1, quad1, timeStep);
  else
    quadtree_insert_lines(tree->quad1, quad1, timeStep);
  }
  if (quad2->head) {
    tree->quad2 = quad_tree_new(xmid, xmax, ymin, ymid);
  if (quad2->num_lines > INSERT_COARSE_LIM)
    cilk_spawn quadtree_insert_lines(tree->quad2, quad2, timeStep);
  else
    quadtree_insert_lines(tree->quad2, quad2, timeStep);
  }
  if (quad3->head) {
    tree->quad3 = quad_tree_new(xmin, xmid, ymid, ymax);
    if (quad3->num_lines > INSERT_COARSE_LIM)
      cilk_spawn quadtree_insert_lines(tree->quad3, quad3, timeStep);
    else
      quadtree_insert_lines(tree->quad3, quad3, timeStep);
  }*/

  if (quad1) {
    tree->quad1 = quad_tree_new(xmin, xmid, ymin, ymid);
    quadtree_insert_lines(tree->quad1, quad1, timeStep, num_quad1);
  }
  if (quad2) {
    tree->quad2 = quad_tree_new(xmid, xmax, ymin, ymid);
    quadtree_insert_lines(tree->quad2, quad2, timeStep, num_quad2);
  }
  if (quad3) {
    tree->quad3 = quad_tree_new(xmin, xmid, ymid, ymax);
    quadtree_insert_lines(tree->quad3, quad3, timeStep, num_quad3);
  }
  if (quad4) {
    tree->quad4 = quad_tree_new(xmid, xmax, ymid, ymax);
    quadtree_insert_lines(tree->quad4, quad4, timeStep, num_quad4);
  }
}
