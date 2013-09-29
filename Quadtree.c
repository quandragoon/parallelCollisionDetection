enum quad_type {Q1 = 1, Q2 = 2, Q3 = 3, Q4 = 4, MUL = 0};

struct line_node {
  struct line_node* next;
  Line line;
};

line_node* line_node_new(Line* line){
  line_node *const new_line = (line_node *const)malloc(sizeof(line_node));
  new_line->line = *line;
  return line_node;  
}

struct line_list {
  size_t num_lines;
  struct line_node* head;
  struct line_node* tail;
};

struct quad_tree {
  struct quad_tree* quad1, quad2, quad3, quad4;
  line_list* lines;
  size_t num_lines; //total lines contained, not the length of 'lines'.
};

typedef struct quad_tree quad_tree;
typedef struct line_node line_node;
typedef struct line_list line_list;

line_list *line_list_new() {
  line_list *const new_list = (line_list *const)malloc(sizeof(line_list));
  new_list->head = NULL;
  new_list->tail = NULL;
  new_list->num_lines = 0;
  return new_list;
}

quad_tree *quad_tree_new() {
  quad_tree *const root = (quad_tree *const)malloc(sizeof(quad_tree));
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
  lines->size++;
  lines->tail->next = NULL;
}



void quadtree_insert_lines(quad_tree* tree, line_list* new_lines) {
  tree->num_lines = new_lines->num_lines;

  if (new_lines->size <= N)
    return;

  line_list* quad1  = line_list_new();
  line_list* quad2  = line_list_new();
  line_list* quad3  = line_list_new();
  line_list* quad4  = line_list_new();
  line_list* parent = line_list_new();

  line_node* cur = new_lines->head;
  quad_type t;
  while (cur->next != NULL){
    t = get_quad_type(tree, cur->line);
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
	return 0;
    }
  }
  tree->quad1 = quad_tree_new();
  tree->quad2 = quad_tree_new();
  tree->quad3 = quad_tree_new();
  tree->quad4 = quad_tree_new();
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
