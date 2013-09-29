struct line_node {
  struct line_node* next;
  Line line;
};

struct line_list {
  struct line_node* head;
  struct line_node* tail;
};

struct quad_tree {
  struct quad_tree* quad1, quad2, quad3, quad4;
  line_list* lines;
};

typedef struct quad_tree quad_tree;
typedef struct line_node line_node;
typedef struct line_list line_list;

line_list *line_list_new() {
  line_list *const new_list = malloc(sizeof(line_list));
  new_list->head = NULL;
  new_list->tail = NULL;
  return new_list;
}

quad_tree *quad_tree_new() {
  quad_tree *const root = malloc(sizeof(quad_tree));
  root->quad1 = root->quad2 = root->quad3 = root->quad4 = NULL;
  root->lines = line_list_new();
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
  lines->tail->next = NULL;
}

void quadtree_insert(quad_tree* tree, line_node* new_line) {
  insert_line(tree->lines, new_line);
}
