//
// Created by armando-albornoz on 2/15/25.
//

#pragma once

#include <cstdint>
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <random>
#include <chrono>
#include <utility>
#include <bits/random.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::vector<Point_2> Points;
typedef std::vector<std::pair<Point_2, Point_2>> Edges;

// Input: arbitrary list of two dimensional points

// Output: Since the output is a polygon, and we can represent a polygon using its vertices,
// the output is a list of vertices of the polygon in clockwise order

// Assume the n points are distinct. First, I am going to implement the SlowConvexHullAlgorithm. To implement this algorithm,
// we can take advantage of primitive operations. Notice that implementing this operations is not neccesary since CGAL provides them


// Slow Convex Hull Algorithm
Points slowConvexHull (const Points& pointSet)
{

  Points convexHull;  // Let convexHull contain the list of points of the convex hull
  Edges E; // The list of edges we have tested so far

 // For each pair of different points
  for (const auto& p : pointSet) {
    for (const auto& q : pointSet) {
      if (p == q) continue;

      bool valid = true;

      // For all points different from p and q perform sidedness test
      for (const auto& r : pointSet) {
        if (r == p || r == q) continue;

        if (CGAL::left_turn(p,q,r)) {
          valid = false;
          break;
        }
      }

      // If the edge is valid add it to E
      if (valid) E.emplace_back(p,q);

    }
  }

  // Construct list of vertices from the valid edges in clockwise order.
  // This is one simple way to do this, but I saw other ways that make use of angle and centroids
  if (!E.empty()) {
    convexHull.push_back(E.at(0).first);
    convexHull.push_back(E.at(0).second);

    for (int i = 1; i < E.size(); i++) {
      if (E[i].first == convexHull.back()) {
        convexHull.push_back(E[i].second);
      }
    }
  }

  return convexHull;
}


Points incrementalAlgorithm (Points& pointSet) {
  // We first sort the point set lexicographically

  std::sort(pointSet.begin(), pointSet.end(),[](const Point_2& p, const Point_2& q) {
    return p.x() < q.x() || (p.x() == q.x() && p.y() < q.y());
  });

  Points upperHull;

  //Put the points p1 and p2 in the upperhull
  upperHull.push_back(pointSet[0]);
  upperHull.push_back(pointSet[1]);

  // Go through sorted points from left to right to build the upper hull
  for (int i= 2; i < pointSet.size(); i++) {
    Point_2 pi = pointSet[i];

    while (upperHull.size() >= 2) {
      // getting the last 3 points (including pi)
        Point_2 a = upperHull[upperHull.size() - 1];
        Point_2 b = upperHull[upperHull.size() - 2];

      // check for right turn
      if (CGAL::left_turn(a,b,pi)) {
        upperHull.pop_back();
      }
      else {
        break;
      }
    }
    upperHull.push_back(pi);
  }

  Points lowerHull;

  lowerHull.push_back(pointSet[pointSet.size() - 1]);
  lowerHull.push_back(pointSet[pointSet.size() - 2]); // at least two points is an assumption

  // Go through sorted points from right to left to build the lower hull
  for (int i = pointSet.size() - 3; i >= 0 ; i--) {
    Point_2 pi = pointSet[i];

    while (lowerHull.size() >= 2) {
      // getting the last 3 points (including pi)
      Point_2 a = lowerHull[lowerHull.size() - 1];
      Point_2 b = lowerHull[lowerHull.size() - 2];

      // check for right turn
      if (CGAL::left_turn(a,b,pi)) {
        lowerHull.pop_back();
      }
      else {
        break;
      }
    }
    lowerHull.push_back(pi);
  }


 // Combine the upper and lower hulls removing duplicates
  Points convexHull;

  convexHull.insert(convexHull.end(), upperHull.begin(), upperHull.end() - 1);
  convexHull.insert(convexHull.end(), lowerHull.begin(), lowerHull.end() - 1);

  return convexHull;
}

// gift wrapping algorithm

Points jarvisMarch(const Points& pointSet) {

  Points convexHull;
  // Find the leftmost point
  Point_2 leftmostPoint = pointSet[0];
  for (const auto& p : pointSet) {
    if(p.x() < leftmostPoint.x()) leftmostPoint = p;
  }

  Point_2 hullPoint = leftmostPoint;
  Point_2 endpoint;

  do {
    convexHull.push_back(hullPoint);
    endpoint = pointSet[0];

     // look for next point on the convexhull
    for (const auto& p : pointSet) {
      if (p == hullPoint) continue;

      // Use primitive sideness operation to find the leftmost point
      if (endpoint == hullPoint || CGAL::left_turn(hullPoint, endpoint, p)) {
        endpoint = p;
      }
    }

    hullPoint = endpoint;
  } while (endpoint != convexHull[0]);

  return convexHull;
}