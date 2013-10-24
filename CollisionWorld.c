/** 
 * CollisionWorld.c -- detect and handle line segment intersections
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

#include "./CollisionWorld.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>

#include "./IntersectionDetection.h"
#include "./IntersectionEventList.h"
#include "./Line.h"
#include "./Quadtree.h"

// Coarsening for computing intersection in parallel
#define INTERSECT_COARSE_LIM 20


CollisionWorld* CollisionWorld_new(const unsigned int capacity) {
  assert(capacity > 0);

  CollisionWorld* collisionWorld = malloc(sizeof(CollisionWorld));
  if (collisionWorld == NULL) {
    return NULL;
  }

  collisionWorld->numLineWallCollisions = 0;
  collisionWorld->numLineLineCollisions = 0;
  collisionWorld->timeStep = 0.5;
  collisionWorld->lines = malloc(capacity * sizeof(Line*));
  collisionWorld->line_nodes = malloc(capacity * sizeof(line_node*));
  collisionWorld->numOfLines = 0;
  return collisionWorld;
}

void CollisionWorld_delete(CollisionWorld* collisionWorld) {
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    free(collisionWorld->lines[i]);
    free(collisionWorld->line_nodes[i]);
  }
  free(collisionWorld->lines);
  /*
  line_node * cur, * prev;
  cur = *collisionWorld->line_nodes;
  while (cur->next){
	  prev = cur;
	  cur = cur->next;
	  free(prev);
  }
  */
  free(collisionWorld->line_nodes);
  free(collisionWorld);
}

unsigned int CollisionWorld_getNumOfLines(CollisionWorld* collisionWorld) {
  return collisionWorld->numOfLines;
}

void CollisionWorld_addLine(CollisionWorld* collisionWorld, Line *line) {
  collisionWorld->lines[collisionWorld->numOfLines] = line;
  collisionWorld->line_nodes[collisionWorld->numOfLines] = line_node_new(line);
  collisionWorld->numOfLines++;
}

Line* CollisionWorld_getLine(CollisionWorld* collisionWorld,
                             const unsigned int index) {
  if (index >= collisionWorld->numOfLines) {
    return NULL;
  }
  return collisionWorld->lines[index];
}

void CollisionWorld_updateLines(CollisionWorld* collisionWorld) {
  CollisionWorld_detectIntersection(collisionWorld);
  CollisionWorld_updatePosition(collisionWorld);
  CollisionWorld_lineWallCollision(collisionWorld);
}

void CollisionWorld_updatePosition(CollisionWorld* collisionWorld) {
  double t = collisionWorld->timeStep;
  Vec displacement;
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *line = collisionWorld->lines[i];
    displacement = Vec_multiply(line->velocity, t);
    line->p1 = Vec_add(line->p1, displacement);
    line->p2 = Vec_add(line->p2, displacement);
  }
}

// Puts all points in the given collision_world into a quad_tree and
// returns the quad_tree.
quad_tree* build_quadtree(CollisionWorld* collision_world) {
  quad_tree* tree = quad_tree_new(BOX_XMIN, BOX_XMAX, BOX_YMIN, BOX_YMAX);
  tree->num_lines = collision_world->numOfLines;

  // Insert all the lines into the root of the tree if total number of
  // lines is less than N
  if (tree->num_lines <= N) {
    for (int i = 0; i < collision_world->numOfLines; ++i) {
      insert_line(&tree->lines, collision_world->line_nodes[i]);
    }
    return tree;
  }

  // quad1, quad2, quad3, quad4 store all line segments that can be completely
  // inserted into smaller quad_trees contained in the given quad_tree.
  //
  // lines contains all line segments that cannot be completely inserted into
  // a single sub-quad_tree of the given quad_tree.
  line_node *quad1, *quad2, *quad3, *quad4, *lines;
  quad1 = quad2 = quad3 = quad4 = lines = NULL;
  int num_quad1, num_quad2, num_quad3, num_quad4;
  num_quad1 = num_quad2 = num_quad3 = num_quad4 = 0;

  // Iterate through all line segments contained in the current quad_tree, and determine
  // which sub-quad_tree a line segment can be inserted into, if any exists.
  int type;
  for (int i = 0; i < collision_world->numOfLines; i++) {
    update_box(collision_world->lines[i], collision_world->timeStep);
    line_node* ptr_node = collision_world->line_nodes[i];
    type = get_quad_type(tree, ptr_node, collision_world->timeStep);
    switch (type) {
      case Q1_TYPE:
        insert_line(&quad1, ptr_node);
        num_quad1++;
        break;
      case Q2_TYPE:
        insert_line(&quad2, ptr_node);
        num_quad2++;
        break;
      case Q3_TYPE:
        insert_line(&quad3, ptr_node);
        num_quad3++;
        break;
      case Q4_TYPE:
        insert_line(&quad4, ptr_node);
        num_quad4++;
        break;
      case MUL_TYPE:
        insert_line(&lines, ptr_node);
        break;
      default:
        return NULL;
    }
  }

  // This is used to determine the dimensions of the sub-quad_trees
  double X_MID = (BOX_XMAX + BOX_XMIN) / 2.0;
  double Y_MID = (BOX_YMAX + BOX_YMIN) / 2.0;

  tree->lines = lines;

  if (quad1) {
    tree->quad1 = quad_tree_new(BOX_XMIN, X_MID, BOX_YMIN, Y_MID);
    cilk_spawn quadtree_insert_lines(tree->quad1, quad1, collision_world->timeStep, num_quad1);
  }
  if (quad2) {
    tree->quad2 = quad_tree_new(X_MID, BOX_XMAX, BOX_YMIN, Y_MID);
    cilk_spawn quadtree_insert_lines(tree->quad2, quad2, collision_world->timeStep, num_quad2);
  }
  if (quad3) {
    tree->quad3 = quad_tree_new(BOX_XMIN, X_MID, Y_MID, BOX_YMAX);
    cilk_spawn quadtree_insert_lines(tree->quad3, quad3, collision_world->timeStep, num_quad3);
  }
  if (quad4) {
    // Note: We do not cilk_spawn here since if we cilk_spawn, the current thread
    // will be wasted waiting for all the other threads to complete
    tree->quad4 = quad_tree_new(X_MID, BOX_XMAX, Y_MID, BOX_YMAX);
    quadtree_insert_lines(tree->quad4, quad4, collision_world->timeStep, num_quad4);
  }

  return tree;
}

// Method that computes the list of intersections within the given quad_tree
IntersectionEventList CollisionWorld_getIntersectionEvents(quad_tree* tree,
  double timeStep, line_node* upstream_lines) {
  IntersectionEventList intersectionEventList = IntersectionEventList_make();
  if (tree == NULL) return intersectionEventList;

  line_node* first_node;
  line_node* second_node;
  first_node = tree->lines;

  // First iterate through all pairs of line segments that cannot be
  // completely inserted into sub-quad_trees
  while (first_node != NULL) {
    second_node = first_node->next;
    while (second_node != NULL) {
      Line* l1 = first_node->line;
      Line* l2 = second_node->line;

      // intersect expects compareLines(l1, l2) < 0 to be true.
      // Swap l1 and l2, if necessary.
      if (compareLines(l1, l2) >= 0) {
        Line *temp = l1;
        l1 = l2;
        l2 = temp;
      }
      IntersectionType intersectionType = intersect(l1, l2, timeStep);
      if (intersectionType != NO_INTERSECTION) {
        IntersectionEventList_appendNode(&intersectionEventList, l1, l2,
                                         intersectionType);
      }
      second_node = second_node->next;
    }
    first_node = first_node->next;
  }

  // Now iterate through all pairs of line segments (a,b) where a is a line
  // segment that cannot be completely inserted into any one sub-quad_tree,
  // and b is a line segment of a quad_tree that contains the current quad_tree
  // as a sub-quad_tree
  first_node = tree->lines;
  while (first_node != NULL) {
    second_node = upstream_lines;
    while (second_node != NULL) {
      Line* l1 = first_node->line;
      Line* l2 = second_node->line;

      // intersect expects compareLines(l1, l2) < 0 to be true.
      // Swap l1 and l2, if necessary.
      if (compareLines(l1, l2) >= 0) {
        Line *temp = l1;
        l1 = l2;
        l2 = temp;
      }

      IntersectionType intersectionType = intersect(l1, l2, timeStep);
      if (intersectionType != NO_INTERSECTION) {
        IntersectionEventList_appendNode(&intersectionEventList, l1, l2,
                                         intersectionType);
      }
      second_node = second_node->next;
    }
    first_node = first_node->next;
  }

  IntersectionEventList intersectionEventListQuad1;
  IntersectionEventList intersectionEventListQuad2;
  IntersectionEventList intersectionEventListQuad3;
  IntersectionEventList intersectionEventListQuad4;

  // We can now propagate tree->lines to all lower sub-quad_trees
  // This does not lead to data races since no element of upstream_lines
  // is modified by the functions that use tree->lines, and tree->lines
  // is a variable that is only being read by multiple threads
  merge_lists(&tree->lines, upstream_lines);

  // For large quad_trees we perform the operation of computing intersections in parallel
  if (tree->num_lines > INTERSECT_COARSE_LIM) {
    intersectionEventListQuad1 = \
      cilk_spawn CollisionWorld_getIntersectionEvents(tree->quad1, timeStep, tree->lines);
    intersectionEventListQuad2 = \
      cilk_spawn CollisionWorld_getIntersectionEvents(tree->quad2, timeStep, tree->lines);
    intersectionEventListQuad3 = \
      cilk_spawn CollisionWorld_getIntersectionEvents(tree->quad3, timeStep, tree->lines);
    intersectionEventListQuad4 = \
      CollisionWorld_getIntersectionEvents(tree->quad4, timeStep, tree->lines);
    cilk_sync;
  } else {
    // For very small quad_trees we do not pay the overhead of spawning
    // new threads
    intersectionEventListQuad1 = \
      CollisionWorld_getIntersectionEvents(tree->quad1, timeStep, tree->lines);
    intersectionEventListQuad2 = \
      CollisionWorld_getIntersectionEvents(tree->quad2, timeStep, tree->lines);
    intersectionEventListQuad3 = \
      CollisionWorld_getIntersectionEvents(tree->quad3, timeStep, tree->lines);
    intersectionEventListQuad4 = \
      CollisionWorld_getIntersectionEvents(tree->quad4, timeStep, tree->lines);
  }

  // Merge the intersections obtained in sub-quad_trees so that we
  // now have one unified list.
  IntersectionEventList_mergeLists(&intersectionEventList,
    &intersectionEventListQuad1);
  IntersectionEventList_mergeLists(&intersectionEventList,
    &intersectionEventListQuad2);
  IntersectionEventList_mergeLists(&intersectionEventList,
    &intersectionEventListQuad3);
  IntersectionEventList_mergeLists(&intersectionEventList,
    &intersectionEventListQuad4);

  return intersectionEventList;
}

void CollisionWorld_detectIntersection(CollisionWorld* collisionWorld) {
  quad_tree* tree = build_quadtree(collisionWorld);
  // Use the constructed quad_tree to detect line-line collisions
  // All line-line intersections are recorded in intersectionEventList
  IntersectionEventList intersectionEventList = \
    CollisionWorld_getIntersectionEvents(tree, collisionWorld->timeStep, NULL);
  collisionWorld->numLineLineCollisions += intersectionEventList.numIntersections;
  quad_tree_delete(tree);
  // Sort the intersection event list.
  IntersectionEventNode* startNode = intersectionEventList.head;
  while (startNode != NULL) {
    IntersectionEventNode* minNode = startNode;
    IntersectionEventNode* curNode = startNode->next;
    while (curNode != NULL) {
      if (IntersectionEventNode_compareData(curNode, minNode) < 0) {
        minNode = curNode;
      }
      curNode = curNode->next;
    }
    if (minNode != startNode) {
      IntersectionEventNode_swapData(minNode, startNode);
    }
    startNode = startNode->next;
  }

  // Call the collision solver for each intersection event.
  IntersectionEventNode* curNode = intersectionEventList.head;

  while (curNode != NULL) {
    CollisionWorld_collisionSolver(collisionWorld, curNode->l1, curNode->l2,
                                   curNode->intersectionType);
    curNode = curNode->next;
  }
  IntersectionEventList_deleteNodes(&intersectionEventList);
}

unsigned int CollisionWorld_getNumLineWallCollisions(
    CollisionWorld* collisionWorld) {
  return collisionWorld->numLineWallCollisions;
}

unsigned int CollisionWorld_getNumLineLineCollisions(
    CollisionWorld* collisionWorld) {
  return collisionWorld->numLineLineCollisions;
}

void CollisionWorld_collisionSolver(CollisionWorld* collisionWorld, Line *l1,
                                    Line *l2, IntersectionType intersectionType) {
  assert(compareLines(l1, l2) < 0);
  assert(
      intersectionType == L1_WITH_L2 || intersectionType == L2_WITH_L1
          || intersectionType == ALREADY_INTERSECTED);

  // Despite our efforts to determine whether lines will intersect ahead
  // of time (and to modify their velocities appropriately), our
  // simplified model can sometimes cause lines to intersect.  In such a
  // case, we compute velocities so that the two lines can get unstuck in
  // the fastest possible way, while still conserving momentum and kinetic
  // energy.
  if (intersectionType == ALREADY_INTERSECTED) {
    Vec p = getIntersectionPoint(l1->p1, l1->p2, l2->p1, l2->p2);

    if (Vec_length(Vec_subtract(l1->p1, p))
        < Vec_length(Vec_subtract(l1->p2, p))) {
      l1->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l1->p2, p)),
                                  Vec_length(l1->velocity));
    } else {
      l1->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l1->p1, p)),
                                  Vec_length(l1->velocity));
    }
    if (Vec_length(Vec_subtract(l2->p1, p))
        < Vec_length(Vec_subtract(l2->p2, p))) {
      l2->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l2->p2, p)),
                                  Vec_length(l2->velocity));
    } else {
      l2->velocity = Vec_multiply(Vec_normalize(Vec_subtract(l2->p1, p)),
                                  Vec_length(l2->velocity));
    }
    return;
  }

  // Compute the collision face/normal vectors.
  Vec face;
  Vec normal;
  if (intersectionType == L1_WITH_L2) {
    Vec v = Vec_makeFromLine(*l2);
    face = Vec_normalize(v);
  } else {
    Vec v = Vec_makeFromLine(*l1);
    face = Vec_normalize(v);
  }
  normal = Vec_orthogonal(face);

  // Obtain each line's velocity components with respect to the collision
  // face/normal vectors.
  double v1Face = Vec_dotProduct(l1->velocity, face);
  double v2Face = Vec_dotProduct(l2->velocity, face);
  double v1Normal = Vec_dotProduct(l1->velocity, normal);
  double v2Normal = Vec_dotProduct(l2->velocity, normal);

  // Compute the mass of each line (we simply use its length).
  double m1 = Vec_length(Vec_subtract(l1->p1, l1->p2));
  double m2 = Vec_length(Vec_subtract(l2->p1, l2->p2));

  // Perform the collision calculation (computes the new velocities along
  // the direction normal to the collision face such that momentum and
  // kinetic energy are conserved).
  double newV1Normal = ((m1 - m2) / (m1 + m2)) * v1Normal
      + (2 * m2 / (m1 + m2)) * v2Normal;
  double newV2Normal = (2 * m1 / (m1 + m2)) * v1Normal
      + ((m2 - m1) / (m2 + m1)) * v2Normal;

  // Combine the resulting velocities.
  l1->velocity = Vec_add(Vec_multiply(normal, newV1Normal),
                         Vec_multiply(face, v1Face));
  l2->velocity = Vec_add(Vec_multiply(normal, newV2Normal),
                         Vec_multiply(face, v2Face));

  return;
}
