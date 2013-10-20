/**
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

#include "./IntersectionDetection.h"

#include <assert.h>
#include <math.h>

#include "./Line.h"
#include "./Vec.h"

#define THRESHOLD 0.001

static inline double max(double x, double y) {
  return x > y ? x : y;
}
static inline double min(double x, double y) {
  return x < y ? x : y;
}

static inline double dabs(double a) {
  if (a < 0)
    return -a;
  return a;
}

// Quick detect if two lines intersect using rectangles
/*
static inline bool rectangles_not_overlap(Vec p11, Vec p12, Vec p21, Vec p22) {
  return (dabs(p11.x - p21.x) > (dabs(p11.x - p12.x) + dabs(p21.x - p22.x)) && \
          dabs(p11.y - p21.y) > (dabs(p11.y - p12.y) + dabs(p21.y - p22.y)));
}
*/

/*
static inline bool rectangles_overlap(Line* l1, Line* l2) {
  return (*l1->min_x < *l2->max_x + THRESHOLD) && (*l1->max_x + THRESHOLD > *l2->min_x) && (*l1->min_y < *l2->max_y + THRESHOLD) && (*l1->max_y + THRESHOLD > *l2->min_y);
}
*/

static inline bool rectangles_overlap(Line* l1, Line* l2) {
  return (l1->l_x < l2->u_x) && (l1->u_x > l2->l_x) && (l1->l_y < l2->u_y) && (l1->u_y > l2->l_y);
}

// Detect if lines l1 and l2 will intersect between now and the next time step.
IntersectionType intersect(Line *l1, Line *l2, double time) {
  assert(compareLines(l1, l2) < 0);

  /*
  if (rectangles_not_overlap(l1->p1, l1->p2, l2->p1, l2->p2))
    return NO_INTERSECTION;
    */

  if (!rectangles_overlap(l1, l2))
    return NO_INTERSECTION;

  Vec velocity;
  Vec p1;
  Vec p2;
  Vec v1 = Vec_makeFromLine(*l1);
  Vec v2 = Vec_makeFromLine(*l2);

  // Get relative velocity.
  velocity = Vec_subtract(l2->velocity, l1->velocity);

  // Get the parallelogram.
  p1 = Vec_add(l2->p1, Vec_multiply(velocity, time));
  p2 = Vec_add(l2->p2, Vec_multiply(velocity, time));


  int num_line_intersections = 0;
  bool top_intersected = false;
  bool bottom_intersected = false;

  if (intersectLines(l1->p1, l1->p2, l2->p1, l2->p2)) {
    return ALREADY_INTERSECTED;
  }
  if (intersectLines(l1->p1, l1->p2, p1, p2)) {
    num_line_intersections++;
  }
  if (intersectLines(l1->p1, l1->p2, p1, l2->p1)) {
    num_line_intersections++;
    top_intersected = true;
  }
  if (intersectLines(l1->p1, l1->p2, p2, l2->p2)) {
    num_line_intersections++;
    bottom_intersected = true;
  }

  if (num_line_intersections == 2) {
    return L2_WITH_L1;
  }

  if (pointInParallelogram(l1->p1, l2->p1, l2->p2, p1, p2)
      && pointInParallelogram(l1->p2, l2->p1, l2->p2, p1, p2)) {
    return L1_WITH_L2;
  }

  if (num_line_intersections == 0) {
    return NO_INTERSECTION;
  }

  double angle = Vec_angle(v1, v2);

  if (top_intersected) {
    if (angle < 0) {
      return L2_WITH_L1;
    } else {
      return L1_WITH_L2;
    }
  }

  if (bottom_intersected) {
    if (angle > 0) {
      return L2_WITH_L1;
    } else {
      return L1_WITH_L2;
    }
  }

  return L1_WITH_L2;
}



// Obtain the intersection point for two intersecting line segments.
Vec getIntersectionPoint(Vec p1, Vec p2, Vec p3, Vec p4) {
  double u;

  u = ((p4.x - p3.x) * (p1.y - p3.y) - (p4.y - p3.y) * (p1.x - p3.x))
      / ((p4.y - p3.y) * (p2.x - p1.x) - (p4.x - p3.x) * (p2.y - p1.y));

  return Vec_add(p1, Vec_multiply(Vec_subtract(p2, p1), u));
}


