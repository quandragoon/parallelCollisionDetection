#include "./Quadtree.h"
#include "./Line.h"
#include "./Vec.h"

#include "assert.h"

line_node* line_node_new(Line* line) {
  line_node *const new_line = (line_node *const)malloc(sizeof(struct line_node));
  new_line->line = line;
  return new_line;
}

line_list *line_list_new() {
  line_list * new_list = (line_list *)malloc(sizeof(struct line_list));
  new_list->head = NULL;
  new_list->tail = NULL;
  new_list->num_lines = 0;
  return new_list;
}

quad_tree *quad_tree_new(double xmin, double xmax, double ymin, double ymax) {
  quad_tree * root = (quad_tree *)malloc(sizeof(quad_tree));
  root->quad1 = root->quad2 = root->quad3 = root->quad4 = NULL;
  root->lines = line_list_new();
  root->num_lines = 0;
  root->xmin = xmin;
  root->xmax = xmax;
  root->ymin = ymin;
  root->ymax = ymax;
  return root;
}

void line_list_delete(line_list* list){
  if (list->head == NULL)
    return;

  line_node *prev, *cur;
  prev = cur = list->head;
  while(cur != NULL){
    cur = cur->next;
    free(prev);
    prev = cur;
  }
}

void quad_tree_delete(quad_tree * tree) {
  //Call recursively
  if (tree->quad1 != NULL)
    quad_tree_delete(tree->quad1);
  if (tree->quad2 != NULL)
    quad_tree_delete(tree->quad2);
  if (tree->quad3 != NULL)
    quad_tree_delete(tree->quad3);
  if (tree->quad4 != NULL)
    quad_tree_delete(tree->quad4);
  line_list_delete(tree->lines);
  free(tree->lines);
}

// Inserts a new line into the given linked list, making sure that
// the input line is not modified by this operation in any way
void insert_line(line_list* lines, line_node* new_line) {
  line_node* node = line_node_new(new_line->line);
  if (lines->head == NULL) {
    lines->head = node;
    lines->tail = node;
  } else {
    lines->tail->next = node;
    lines->tail = node;
  }
  lines->num_lines++;
  lines->tail->next = NULL;
}

// Merges the two lists, does not modify list2
void merge_lists(line_list* list1, line_list* list2) {
  if (list2 == NULL) return;
  if (list1->head == NULL) {
    list1->head = list2->head;
    list1->tail = list2->tail;
  } else {
    list1->tail->next = list2->head;
    if (list2->tail != NULL)
      list1->tail = list2->tail;
  }
  list1->num_lines += list2->num_lines;
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

void quadtree_insert_lines(quad_tree* tree, line_list* new_lines, double timeStep) {
  tree->num_lines = new_lines->num_lines;
  double xmax = tree->xmax;
  double xmin = tree->xmin;
  double ymax = tree->ymax;
  double ymin = tree->ymin;

  if (new_lines->num_lines <= N) {
    tree->lines = new_lines;
    return;
  }

  line_list* quad1  = line_list_new();
  line_list* quad2  = line_list_new();
  line_list* quad3  = line_list_new();
  line_list* quad4  = line_list_new();
  line_list* parent = line_list_new();

  line_node* cur = new_lines->head;
  int type;
  while (cur != NULL) {
    type = get_quad_type(tree, cur, timeStep);
    switch (type) {
      case Q1_TYPE:
        insert_line(quad1, cur);
        break;
      case Q2_TYPE:
        insert_line(quad2, cur);
        break;
      case Q3_TYPE:
        insert_line(quad3, cur);
        break;
      case Q4_TYPE:
        insert_line(quad4, cur);
        break;
      case MUL_TYPE:
        insert_line(parent, cur);
        break;
      default:
        return;
    }
    cur = cur->next;
  }

  assert(new_lines->num_lines == \
    (parent->num_lines + quad1->num_lines + quad2->num_lines \
      + quad3->num_lines + quad4->num_lines));

  double xmid = (xmin + xmax) / 2.0;
  double ymid = (ymin + ymax) / 2.0;
  tree->lines = parent;
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

  if (quad1->head) {
    tree->quad1 = quad_tree_new(xmin, xmid, ymin, ymid);
    quadtree_insert_lines(tree->quad1, quad1, timeStep);
  }
  if (quad2->head) {
    tree->quad2 = quad_tree_new(xmid, xmax, ymin, ymid);
    quadtree_insert_lines(tree->quad2, quad2, timeStep);
  }
  if (quad3->head) {
    tree->quad3 = quad_tree_new(xmin, xmid, ymid, ymax);
    quadtree_insert_lines(tree->quad3, quad3, timeStep);
  }
  if (quad4->head) {
    tree->quad4 = quad_tree_new(xmid, xmax, ymid, ymax);
    quadtree_insert_lines(tree->quad4, quad4, timeStep);
  }
}