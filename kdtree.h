    #include <CGAL/Simple_cartesian.h>
    #include <vector>
    #include <algorithm>
    #include <iostream>
    #include <limits>

    typedef CGAL::Simple_cartesian<double> Kernel;
    typedef Kernel::Point_2 Point_2;
    typedef Kernel::Iso_rectangle_2 Iso_rectangle_2;

    struct KDTreeNode {
        Point_2 point;
        double split_value;
        KDTreeNode* left;
        KDTreeNode* right;
        KDTreeNode(Point_2 p, double split = 0) : point(p), split_value(split), left(nullptr), right(nullptr) {}
    };

    struct Region {
        double xmin, xmax;
        double ymin, ymax;

        Region(double xmin = -std::numeric_limits<double>::infinity(),
               double xmax = std::numeric_limits<double>::infinity(),
               double ymin = -std::numeric_limits<double>::infinity(),
               double ymax = std::numeric_limits<double>::infinity())
            : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax) {}

        Iso_rectangle_2 to_iso_rectangle() const {
            return Iso_rectangle_2(Point_2(xmin, ymin), Point_2(xmax, ymax));
        }
    };

    KDTreeNode* buildKDTree(std::vector<Point_2>& points, int depth = 0) {
        if (points.empty()) return nullptr;

        int axis = depth % 2; // 0 for x-axis, 1 for y-axis

        auto cmp = [axis](const Point_2& a, const Point_2& b) {
            return (axis == 0) ? (a.x() < b.x()) : (a.y() < b.y());
        };
        std::sort(points.begin(), points.end(), cmp);

        size_t median = points.size() / 2;
        double split_value = (axis == 0) ? points[median].x() : points[median].y();

        KDTreeNode* node = new KDTreeNode(points[median], split_value);

        std::vector<Point_2> leftPoints(points.begin(), points.begin() + median);
        std::vector<Point_2> rightPoints(points.begin() + median + 1, points.end());

        node->left = buildKDTree(leftPoints, depth + 1);
        node->right = buildKDTree(rightPoints, depth + 1);

        return node;
    }

    void searchKDTree(KDTreeNode* node, const Region& current_region, const Iso_rectangle_2& query_range, std::vector<Point_2>& result, int depth = 0) {
        if (node == nullptr) return;

        Iso_rectangle_2 current_rect = current_region.to_iso_rectangle();
        if (!CGAL::do_intersect(current_rect, query_range)) {
            return;
        }

        // Manual check for point inclusion in the query range
        const Point_2& point = node->point;
        Point_2 q_min = query_range.min();
        Point_2 q_max = query_range.max();
        double x = point.x();
        double y = point.y();

        if (x >= q_min.x() && x <= q_max.x() && y >= q_min.y() && y <= q_max.y()) {
            result.push_back(point);
        }

        int axis = depth % 2;
        Region left_region = current_region;
        Region right_region = current_region;

        if (axis == 0) {
            left_region.xmax = node->split_value;
            right_region.xmin = node->split_value;
        } else {
            left_region.ymax = node->split_value;
            right_region.ymin = node->split_value;
        }

        searchKDTree(node->left, left_region, query_range, result, depth + 1);
        searchKDTree(node->right, right_region, query_range, result, depth + 1);
    }

