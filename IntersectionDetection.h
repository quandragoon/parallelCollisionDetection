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

#ifndef INTERSECTIONDETECTION_H_
#define INTERSECTIONDETECTION_H_

#include "./Line.h"
#include "./Vec.h"

typedef enum {
  NO_INTERSECTION,
  L1_WITH_L2,
  L2_WITH_L1,
  ALREADY_INTERSECTED
} IntersectionType;

// Detect if line l1 and l2 will be intersected in the next time step.
// Precondition: compareLines(l1, l2) < 0 must be true.
IntersectionType intersect(Line *l1, Line *l2, double time);

// Calculate the cross product.
static inline double crossProduct(double x1, double y1, double x2, double y2) {
  return x1 * y2 - x2 * y1;
}

// Check the direction of two lines (pi, pj) and (pi, pk).
static inline double direction(Vec pi, Vec pj, Vec pk) {
  return crossProduct(pk.x - pi.x, pk.y - pi.y, pj.x - pi.x, pj.y - pi.y);
}

// Check if a point is in the parallelogram.
static inline bool pointInParallelogram(Vec point, Vec p1, Vec p2, Vec p3, Vec p4) {
  double d1 = direction(p1, p2, point);
  double d2 = direction(p3, p4, point);
  if ((d1 < 0 && d2 < 0) || (d1 > 0 && d2 > 0)) return false;
  double d3 = direction(p1, p3, point);
  double d4 = direction(p2, p4, point);
  if ((d3 < 0 && d4 < 0) || (d3 > 0 && d4 > 0)) return false;

  return true;
}



// Check if a point pk is in the line segment (pi, pj).
// pi, pj, and pk must be collinear.
static inline bool onSegment(Vec pi, Vec pj, Vec pk) {
  if (((pi.x <= pk.x && pk.x <= pj.x) || (pj.x <= pk.x && pk.x <= pi.x))
      && ((pi.y <= pk.y && pk.y <= pj.y) || (pj.y <= pk.y && pk.y <= pi.y))) {
    return true;
  }
  return false;
}

// Check if a point P is on one side of the line EF or another
static inline bool which_side(Vec E, Vec F, Vec P) {
  return (F.x - E.x) * (P.y - F.y) - (F.y - E.y) * (P.x - F.x) >= 0;
}

// Check if two lines intersect.
static inline bool intersectLines(Vec p1, Vec p2, Vec p3, Vec p4) {
  return which_side(p1, p2, p3) != which_side(p1, p2, p4) && \
    which_side(p3, p4, p1) != which_side(p3, p4, p2);
}

// Obtain the intersection point for two intersecting line segments.
Vec getIntersectionPoint(Vec p1, Vec p2, Vec p3, Vec p4);

#endif  // INTERSECTIONDETECTION_H_
