#include <fstream>
#include <vector>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/convex_hull_2.h>
#include <random>
#include <cmath>
#include <chrono>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Iso_rectangle_2 Iso_rectangle_2;

// Include the range tree and kd-tree implementations
#include "rangetree.h"
#include "kdtree.h"

// Data generation functions (unchanged)
std::vector<Point_2> generateRandomPoints(int n) {
    CGAL::Random rng;
    CGAL::Random_points_in_square_2<Point_2> gen(100.0, rng);
    std::vector<Point_2> points;
    CGAL::cpp11::copy_n(gen, n, std::back_inserter(points));
    return points;
}

std::vector<Point_2> generateConvexPoints(int n) {
    auto points = generateRandomPoints(n);
    std::vector<Point_2> convexHull;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(convexHull));
    return convexHull;
}

std::vector<Point_2> generateCirclePoints(int n, double radius = 50.0) {
    std::vector<Point_2> points;
    double angleIncrement = 2 * M_PI / n;
    for (int i = 0; i < n; ++i) {
        double angle = i * angleIncrement;
        double x = radius * cos(angle) + 50.0;
        double y = radius * sin(angle) + 50.0;
        points.emplace_back(x, y);
    }
    return points;
}

std::vector<Point_2> generateGaussianPoints(int n, double meanX = 50.0, double meanY = 50.0, double stddev = 10.0) {
    std::vector<Point_2> points;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> distX(meanX, stddev);
    std::normal_distribution<> distY(meanY, stddev);

    for (int i = 0; i < n; ++i) {
        double x = distX(gen);
        double y = distY(gen);
        points.emplace_back(x, y);
    }
    return points;
}


// Timing function
template<typename Func>
double measureTime(Func func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count();
}

// Query range function
Iso_rectangle_2 getQueryRange(const std::vector<Point_2>& points) {
    double xmin = 40.0, xmax = 60.0;
    double ymin = 40.0, ymax = 60.0;
    return Iso_rectangle_2(Point_2(xmin, ymin), Point_2(xmax, ymax));
}

// Log results to CSV
void logToCSV(int n, const std::string& distribution, double rangeTreeBuildTime, double rangeTreeQueryTime, double kdTreeBuildTime, double kdTreeQueryTime) {
    std::ofstream csvFile("results.csv", std::ios::app);
    csvFile << n << "," << distribution << ","
            << rangeTreeBuildTime << "," << rangeTreeQueryTime << ","
            << kdTreeBuildTime << "," << kdTreeQueryTime << "\n";
    csvFile.close();
}

// Run experiment for a given size and distribution
void runExperiment(int n, const std::string& distribution) {
    std::vector<Point_2> points;

    if (distribution == "random") {
        points = generateRandomPoints(n);
    } else if (distribution == "convex") {
        points = generateConvexPoints(n);
    } else if (distribution == "circle") {
        points = generateCirclePoints(n);
    } else if (distribution == "gaussian") {
        points = generateGaussianPoints(n);
    }

    if (points.empty()) {
        std::cerr << "Error: Empty points vector for distribution " << distribution << " (n=" << n << ")\n";
        return;
    }

    // Measure range tree performance
    RangeTreeNode* rangeTreeRoot = nullptr;
    auto rangeTreeBuildTime = measureTime([&]() { rangeTreeRoot = build2DRangeTree(points); });

    if (!rangeTreeRoot) {
        std::cerr << "Error: Failed to build range tree for distribution " << distribution << " (n=" << n << ")\n";
        return;
    }

    Iso_rectangle_2 queryRange = getQueryRange(points);
    std::vector<Point_2> rangeTreeResult;
    auto rangeTreeQueryTime = measureTime([&]() {
        query2DRangeTree(
            rangeTreeRoot,
            queryRange.xmin(), queryRange.xmax(),
            queryRange.ymin(), queryRange.ymax(),
            rangeTreeResult
        );
    });

    // Measure kd-tree performance
    KDTreeNode* kdTreeRoot = nullptr;
    auto kdTreeBuildTime = measureTime([&]() { kdTreeRoot = buildKDTree(points); });

    if (!kdTreeRoot) {
        std::cerr << "Error: Failed to build kd-tree for distribution " << distribution << " (n=" << n << ")\n";
        return;
    }

    std::vector<Point_2> kdTreeResult;
    Region initial_region; // Creates region (-∞, ∞) in both dimensions
    auto kdTreeQueryTime = measureTime([&]() {
        searchKDTree(kdTreeRoot, initial_region, queryRange, kdTreeResult, 0);
    });

    // Write results to CSV
    logToCSV(n, distribution, rangeTreeBuildTime, rangeTreeQueryTime, kdTreeBuildTime, kdTreeQueryTime);
}

int main() {
    std::vector<int> sizes = {100, 1000, 10000, 100000};
    std::vector<std::string> distributions = {"random", "convex", "circle", "gaussian"};

    std::ofstream csvFile("results.csv");
    csvFile << "Size,Distribution,RangeTree Build Time (s),RangeTree Query Time (s),KDTree Build Time (s),KDTree Query Time (s)\n";
    csvFile.close();

    for (int n : sizes) {
        for (const std::string& distribution : distributions) {
            runExperiment(n, distribution);
        }
    }

    return 0;
}