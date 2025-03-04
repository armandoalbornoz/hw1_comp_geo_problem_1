#pragma once
#include <string>
#include <algorithm>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/partition_2.h>
#include <stack>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <cassert>
#include <list>

#include <CGAL/Partition_traits_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Partition_traits_2<K> Traits1;
typedef Traits1::Point_2 Point_2;
typedef Traits1::Polygon_2 Polygon_2;
typedef std::vector<Polygon_2> Triang;

typedef CGAL::Arr_segment_traits_2<K> Traits;
typedef CGAL::Arr_default_dcel<Traits> Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel> Arrangement;
typedef Arrangement::Vertex_handle Vertex_handle;
typedef Arrangement::Halfedge_handle Halfedge_handle;
typedef Arrangement::Face_handle Face_handle;

// Helper subroutines

// Check if a vertex is an ear
bool is_ear_1(const Polygon_2& polygon, typename Polygon_2::Vertex_const_iterator current_it) {
    auto prev_it = (current_it == polygon.vertices_begin()) ? std::prev(polygon.vertices_end()) : std::prev(current_it);
    auto next_it = (current_it == std::prev(polygon.vertices_end())) ? polygon.vertices_begin() : std::next(current_it);

    // Make a triangle
    Polygon_2 triangle;
    triangle.push_back(*prev_it);
    triangle.push_back(*current_it);
    triangle.push_back(*next_it);

    for (auto it = polygon.vertices_begin(); it != polygon.vertices_end(); ++it) {
        if (it == prev_it || it == current_it || it == next_it) continue;

        // If any other vertex lies inside the triangle, the triangle can't be an ear
        if (CGAL::bounded_side_2(triangle.vertices_begin(), triangle.vertices_end(), *it) != CGAL::ON_UNBOUNDED_SIDE) {
            return false;
        }
    }
    return true;
}

// Check if a vertex is an ear (optimized for convex vertices)
bool is_ear_2(const Polygon_2& polygon, typename Polygon_2::Vertex_const_iterator current_it, const std::vector<typename Polygon_2::Vertex_const_iterator>& concave_vertices) {
    auto prev_it = (current_it == polygon.vertices_begin()) ? std::prev(polygon.vertices_end()) : std::prev(current_it);
    auto next_it = (current_it == std::prev(polygon.vertices_end())) ? polygon.vertices_begin() : std::next(current_it);

    // Make a triangle
    Polygon_2 triangle;
    triangle.push_back(*prev_it);
    triangle.push_back(*current_it);
    triangle.push_back(*next_it);

    // Check if any concave vertex lies inside the triangle
    for (auto it : concave_vertices) {
        if (it == prev_it || it == current_it || it == next_it) continue;
        if (CGAL::bounded_side_2(triangle.vertices_begin(), triangle.vertices_end(), *it) != CGAL::ON_UNBOUNDED_SIDE) {
            return false;
        }
    }
    return true;
}

// Check if a vertex is convex
bool is_convex(const Polygon_2& polygon, typename Polygon_2::Vertex_const_iterator current_it) {
    auto prev_it = (current_it == polygon.vertices_begin()) ? std::prev(polygon.vertices_end()) : std::prev(current_it);
    auto next_it = (current_it == std::prev(polygon.vertices_end())) ? polygon.vertices_begin() : std::next(current_it);

    // Compute the orientation of the three points
    CGAL::Orientation orient = CGAL::orientation(*prev_it, *current_it, *next_it);

    // Convex if the orientation is counterclockwise
    return (orient == CGAL::LEFT_TURN);
}

// Divide vertices into convex and concave
void divide_vertices(const Polygon_2& polygon, std::vector<typename Polygon_2::Vertex_const_iterator>& convex_vertices, std::vector<typename Polygon_2::Vertex_const_iterator>& concave_vertices) {
    for (auto it = polygon.vertices_begin(); it != polygon.vertices_end(); ++it) {
        if (is_convex(polygon, it)) {
            convex_vertices.push_back(it); // Convex vertex
        } else {
            concave_vertices.push_back(it); // Concave vertex
        }
    }
}

// First Algorithm
Triang triangulate_1(Polygon_2 polygon) {
    Triang triangulation;

    while (polygon.size() > 3) {
        for (auto it = polygon.vertices_begin(); it != polygon.vertices_end(); ++it) {
            if (is_ear_1(polygon, it)) {
                // Add the ear triangle to the triangulation
                Polygon_2 triangle;
                auto prev_it = (it == polygon.vertices_begin()) ? std::prev(polygon.vertices_end()) : std::prev(it);
                auto next_it = (it == std::prev(polygon.vertices_end())) ? polygon.vertices_begin() : std::next(it);
                triangle.push_back(*prev_it);
                triangle.push_back(*it);
                triangle.push_back(*next_it);
                triangulation.push_back(triangle);

                // Remove the ear vertex from the polygon
                polygon.erase(it);
                break;
            }
        }
    }
    // Add the last triangle
    triangulation.push_back(polygon);
    return triangulation;
}

// Second Algorithm
Triang triangulate_2(Polygon_2 polygon) {
    Triang triangulation;
    std::vector<typename Polygon_2::Vertex_const_iterator> convex_vertices, concave_vertices;

    // Divide vertices as convex or concave
    divide_vertices(polygon, convex_vertices, concave_vertices);

    while (polygon.size() > 3) {
        for (auto it : convex_vertices) {
            if (is_ear_2(polygon, it, concave_vertices)) {
                // Add the triangle to the triangulation
                Polygon_2 triangle;
                auto prev_it = (it == polygon.vertices_begin()) ? std::prev(polygon.vertices_end()) : std::prev(it);
                auto next_it = (it == std::prev(polygon.vertices_end())) ? polygon.vertices_begin() : std::next(it);
                triangle.push_back(*prev_it);
                triangle.push_back(*it);
                triangle.push_back(*next_it);
                triangulation.push_back(triangle);

                // Remove the ear vertex from the polygon
                polygon.erase(it);

                // Update the convex and concave lists by re-dividing
                convex_vertices.clear();
                concave_vertices.clear();
                divide_vertices(polygon, convex_vertices, concave_vertices);
                break;
            }
        }
    }

    // Add the last triangle
    triangulation.push_back(polygon);
    return triangulation;
}


void TriangulateMonotonePolygon(Arrangement& D) {
    std::vector<Point_2> vertices;
    for (auto it = D.vertices_begin(); it != D.vertices_end(); ++it) {
        vertices.push_back(it->point());
    }
    std::sort(vertices.begin(), vertices.end(), [](const Point_2& a, const Point_2& b) {
        return a.y() > b.y() || (a.y() == b.y() && a.x() < b.x());
    });
    std::vector<bool> is_left_chain(vertices.size(), false);
    Point_2 top_vertex = vertices[0];
    Point_2 bottom_vertex = vertices.back();
    for (int i = 0; i < vertices.size(); i++) {
        if (vertices[i].y() == top_vertex.y() || vertices[i].y() == bottom_vertex.y()) {
            is_left_chain[i] = true;
        } else if (vertices[i].x() <= top_vertex.x() && vertices[i].x() <= bottom_vertex.x()) {
            is_left_chain[i] = true;
        } else {
            is_left_chain[i] = false;
        }
    }
    std::stack<size_t> S;
    S.push(0);
    S.push(1);
    for (int j = 2; j < vertices.size() - 1; ++j) {
        size_t top_index = S.top();
        bool on_different_chains = (is_left_chain[j] != is_left_chain[top_index]);
        if (on_different_chains) {
            while (!S.empty()) {
                size_t popped_index = S.top();
                S.pop();
                if (!S.empty()) {
                    insert(D, Traits::Segment_2(vertices[j], vertices[popped_index]));
                    insert(D, Traits::Segment_2(vertices[j], vertices[S.top()]));
                }
            }
            S.push(j - 1);
            S.push(j);
        } else {
            size_t popped_index = S.top();
            S.pop();
            while (!S.empty() && CGAL::orientation(vertices[S.top()], vertices[popped_index], vertices[j]) == CGAL::LEFT_TURN) {
                CGAL::insert(D, Traits::Segment_2(vertices[j], vertices[popped_index]));
                CGAL::insert(D, Traits::Segment_2(vertices[j], vertices[S.top()]));
                popped_index = S.top();
                S.pop();
            }
            S.push(popped_index);
            S.push(j);
        }
    }
    Point_2 u_n = vertices.back();
    while (S.size() > 2) {
        size_t popped_index = S.top();
        S.pop();
        insert(D, Traits::Segment_2(u_n, vertices[popped_index]));
        insert(D, Traits::Segment_2(u_n, vertices[S.top()]));
    }
}

void triangulatePolygon(const Polygon_2& polygon) {
    if (CGAL::is_y_monotone_2(polygon.vertices_begin(), polygon.vertices_end())) {
        Arrangement D;
        for (auto edge = polygon.edges_begin(); edge != polygon.edges_end(); ++edge) {
            insert(D, *edge);
        }
        TriangulateMonotonePolygon(D);
    } else {
        std::vector<Polygon_2> monotone_polygons;
        Polygon_2 poly_vector(polygon.vertices_begin(), polygon.vertices_end());
        CGAL::y_monotone_partition_2(poly_vector.vertices_begin(), poly_vector.vertices_end(), std::back_inserter(monotone_polygons));
        for (const auto& monotone_polygon : monotone_polygons) {
            Arrangement D;
            for (auto edge = monotone_polygon.edges_begin(); edge != monotone_polygon.edges_end(); ++edge) {
                insert(D, *edge);
            }
            TriangulateMonotonePolygon(D);
        }
    }
}