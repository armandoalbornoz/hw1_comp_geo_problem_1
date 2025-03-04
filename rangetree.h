#ifndef RANGETREE_H
#define RANGETREE_H

#include <CGAL/Simple_cartesian.h>
#include <vector>
#include <algorithm>
#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Iso_rectangle_2 Iso_rectangle_2;

struct RangeTreeNode {
    double split_value;
    RangeTreeNode* left;
    RangeTreeNode* right;
    std::vector<Point_2> points;

    RangeTreeNode(double split_value = 0)
        : split_value(split_value), left(nullptr), right(nullptr) {}
};

inline RangeTreeNode* buildAssociatedStructure(std::vector<Point_2>& points) {
    if (points.empty()) return nullptr;

    auto compareY = [](const Point_2& a, const Point_2& b) { return a.y() < b.y(); };
    std::sort(points.begin(), points.end(), compareY);

    RangeTreeNode* node = new RangeTreeNode();
    node->points = points;
    return node;
}

inline RangeTreeNode* build2DRangeTree(std::vector<Point_2>& points) {
    if (points.empty()) return nullptr;

    auto compareX = [](const Point_2& a, const Point_2& b) { return a.x() < b.x(); };
    std::sort(points.begin(), points.end(), compareX);

    if (points.size() == 1) {
        RangeTreeNode* node = new RangeTreeNode(points[0].x());
        node->points = points;
        return node;
    }

    size_t median = (points.size() - 1) / 2;
    double x_mid = points[median].x();

    std::vector<Point_2> leftPoints(points.begin(), points.begin() + median + 1);
    std::vector<Point_2> rightPoints(points.begin() + median + 1, points.end());

    RangeTreeNode* leftChild = build2DRangeTree(leftPoints);
    RangeTreeNode* rightChild = build2DRangeTree(rightPoints);

    RangeTreeNode* node = new RangeTreeNode(x_mid);
    node->left = leftChild;
    node->right = rightChild;
    node->points = points;

    return node;
}

inline void queryAssociatedStructure(RangeTreeNode* node, double y1, double y2, std::vector<Point_2>& result) {
    if (node == nullptr) return;

    auto low = std::lower_bound(node->points.begin(), node->points.end(), Point_2(0, y1),
        [](const Point_2& a, const Point_2& b) { return a.y() < b.y(); });
    auto high = std::upper_bound(node->points.begin(), node->points.end(), Point_2(0, y2),
        [](const Point_2& a, const Point_2& b) { return a.y() < b.y(); });

    for (auto it = low; it != high; ++it) {
        result.push_back(*it);
    }
}

inline void collectCanonicalSubsets(RangeTreeNode* node, double x1, double x2, std::vector<RangeTreeNode*>& canonicalNodes) {
    RangeTreeNode* splitNode = node;
    while (splitNode && (splitNode->left || splitNode->right)) {
        if (x2 < splitNode->split_value) {
            splitNode = splitNode->left;
        } else if (x1 > splitNode->split_value) {
            splitNode = splitNode->right;
        } else {
            break;
        }
    }
    if (!splitNode) return;

    RangeTreeNode* v = splitNode->left;
    while (v && (v->left || v->right)) {
        if (x1 <= v->split_value) {
            if (v->right) canonicalNodes.push_back(v->right);
            v = v->left;
        } else {
            v = v->right;
        }
    }
    if (v && v->points.back().x() >= x1 && v->points[0].x() <= x2) {
        canonicalNodes.push_back(v);
    }

    v = splitNode->right;
    while (v && (v->left || v->right)) {
        if (x2 >= v->split_value) {
            if (v->left) canonicalNodes.push_back(v->left);
            v = v->right;
        } else {
            v = v->left;
        }
    }
    if (v && v->points.back().x() >= x1 && v->points[0].x() <= x2) {
        canonicalNodes.push_back(v);
    }
}

inline void query2DRangeTree(RangeTreeNode* root, double x1, double x2, double y1, double y2, std::vector<Point_2>& result) {
    if (!root) return; // Critical null check added

    std::vector<RangeTreeNode*> canonicalNodes;
    collectCanonicalSubsets(root, x1, x2, canonicalNodes);

    for (RangeTreeNode* node : canonicalNodes) {
        queryAssociatedStructure(node, y1, y2, result);
    }

    if (root->left == nullptr && root->right == nullptr) {
        const Point_2& point = root->points[0];
        if (point.x() >= x1 && point.x() <= x2 && point.y() >= y1 && point.y() <= y2) {
            result.push_back(point);
        }
    }
}

#endif // RANGETREE_H