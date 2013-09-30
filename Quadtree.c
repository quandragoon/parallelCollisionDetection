#include "./Quadtree.h"
#include "./Line.h"
#include "./Vec.h"
#include "assert.h"
#include <stdlib.h>

/*struct line_node {
  struct line_node* next;
  Line line;
};*/

line_node* line_node_new(Line* line){
  line_node *const new_line = (line_node *const)malloc(sizeof(struct line_node));
  new_line->line = *line;
  return new_line;  
}

/*struct line_list {
  size_t num_lines;
  struct line_node* head;
  struct line_node* tail;
};

struct quad_tree {
  struct quad_tree* quad1, quad2, quad3, quad4;
  line_list* lines;
  size_t num_lines; //total lines contained, not the length of 'lines'.
  double xmin, xmax, ymin, ymax;
};

typedef struct quad_tree quad_tree;
typedef struct line_node line_node;
typedef struct line_list line_list;
*/
line_list *line_list_new() {
  line_list *const new_list = (line_list *const)malloc(sizeof(struct line_list));
  new_list->head = NULL;
  new_list->tail = NULL;
  new_list->num_lines = 0;
  return new_list;
}

quad_tree *quad_tree_new(double xmin, double xmax, double ymin, double ymax) {
  quad_tree *const root = (quad_tree *const)malloc(sizeof(quad_tree));
  root->quad1 = root->quad2 = root->quad3 = root->quad4 = NULL;
  root->lines = line_list_new();
  root->xmin = xmin;
  root->xmax = xmax;
  root->ymin = ymin;
  root->ymax = ymax;
  return root;
}

void insert_line(line_list* lines, line_node* new_line) {
  if (lines->head == NULL) {
    lines->head = new_line;
    lines->tail = new_line;
  } else {
    lines->tail->next = new_line;
    lines->tail = new_line;
  }
  lines->num_lines++;
  lines->tail->next = NULL;
}

quad_type get_quad_type(quad_tree* tree, line_node* node) {
  Line line = node->line;
  Vec p1 = line.p1;
  Vec p2 = line.p2;

  double xmin = tree->xmin;
  double xmax = tree->xmax;
  double ymin = tree->ymin;
  double ymax = tree->ymax;
  double xmid = (xmin + xmax) / 2.0;
  double ymid = (ymin + ymax) / 2.0;

  assert(p1.x > xmin && p1.x < xmax && p1.y > ymin && p1.y < ymax);
  assert(p2.x > xmin && p2.x < xmax && p2.y > ymin && p2.y < ymax);

  if (!(((p1.x - xmid)*(p2.x - xmid) > 0) && ((p1.y - ymid)*(p2.y - ymid)))) {
    return MUL;
  }
  int xid = (p1.x - xmid > 0) ? 1 : 0;
  int yid = (p1.y - ymid > 0) ? 1 : 0;
  int quad = 2 * yid + xid + 1;
  return (quad_type) quad;
}

void quadtree_insert_lines(quad_tree* tree, line_list* new_lines) {
  tree->num_lines = new_lines->num_lines;
  double xmax = tree->xmax;
  double xmin = tree->xmin;
  double ymax = tree->ymax;
  double ymin = tree->ymin;

  if (new_lines->num_lines <= N)
    return;

  line_list* quad1  = line_list_new();
  line_list* quad2  = line_list_new();
  line_list* quad3  = line_list_new();
  line_list* quad4  = line_list_new();
  line_list* parent = line_list_new();

  line_node* cur = new_lines->head;
  quad_type t;
  while (cur->next != NULL){
    t = get_quad_type(tree, cur);
    switch (t){
      case Q1:
	insert_line(quad1, cur);
	break;
      case Q2:
	insert_line(quad2, cur);
	break;
      case Q3:
	insert_line(quad3, cur);
	break;
      case Q4:
	insert_line(quad4, cur);
	break;
      case MUL:
	insert_line(parent, cur);
	break;
      default:
	return;
    }
  }
  double xmid = (xmin + xmax) / 2.0;
  double ymid = (ymin + ymax) / 2.0;
  tree->quad1 = quad_tree_new(xmin, xmid, ymin, ymid);
  tree->quad2 = quad_tree_new(xmid, xmax, ymin, ymid);
  tree->quad3 = quad_tree_new(xmin, xmid, ymid, ymax);
  tree->quad4 = quad_tree_new(xmid, xmax, ymid, ymax);
  quadtree_insert_line_list(tree, parent);
  quadtree_insert_lines(tree->quad1, quad1);
  quadtree_insert_lines(tree->quad2, quad2);
  quadtree_insert_lines(tree->quad3, quad3);
  quadtree_insert_lines(tree->quad4, quad4);
}

//insert a line into a quadrant that it belongs to
void quadtree_insert_line_list(quad_tree* tree, line_list* new_line) {
  tree->lines = new_line;
}
