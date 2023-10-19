#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;

struct Point {
    int id;
    double x, y;
};

double euclidean_distance(const Point& a, const Point& b) {
    double xd = a.x - b.x;
    double yd = a.y - b.y;
    return round(sqrt(xd * xd + yd * yd));
}

double att_distance(const Point& a, const Point& b) {
    double xd = a.x - b.x;
    double yd = a.y - b.y;
    double rij = sqrt((xd * xd + yd * yd) / 10.0);
    double tij = round(rij);
    if (tij < rij) {
        return tij + 1;
    }
    return tij;
}

vector<int> nearest_insertion(const vector<Point>& points, const char* distance_type) {
    vector<int> tour;
    tour.push_back(points[0].id);

    for (size_t i = 1; i < points.size(); ++i) {
        double min_dist = 1e9;
        int insert_position = -1;

        for (size_t j = 0; j < tour.size(); ++j) {
            double dist;
            const Point& current_point = points[tour[j] - 1];
            const Point& next_point = (j < tour.size() - 1) ? points[tour[j + 1] - 1] : points[tour[0] - 1];

            if (strcmp(distance_type, "EUC_2D") == 0) {
                dist = euclidean_distance(points[i], current_point) + euclidean_distance(points[i], next_point) - euclidean_distance(current_point, next_point);
            } else {
                dist = att_distance(points[i], current_point) + att_distance(points[i], next_point) - att_distance(current_point, next_point);
            }

            if (dist < min_dist) {
                min_dist = dist;
                insert_position = j;
            }
        }

        tour.insert(tour.begin() + insert_position + 1, points[i].id);
    }

    return tour;
}

int main() {
    int dimension;
    char edge_weight_type[100];
    vector<Point> points;

    string line;
    while (getline(cin, line)) {
        if (line.find("DIMENSION") != string::npos) {
            sscanf(line.substr(line.find(":") + 1).c_str(), "%d", &dimension);
        } else if (line.find("EDGE_WEIGHT_TYPE") != string::npos) {
            sscanf(line.substr(line.find(":") + 1).c_str(), "%s", edge_weight_type);
        } else if (line.find("NODE_COORD_SECTION") != string::npos) {
            for (int i = 0; i < dimension; i++) {
                Point point;
                cin >> point.id >> point.x >> point.y;
                points.push_back(point);
            }
        } else if (line.find("EOF") != string::npos) {
            break;
        }
    }

    vector<int> tour = nearest_insertion(points, edge_weight_type);

    double total_distance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        const Point& p1 = points[tour[i] - 1];
        const Point& p2 = points[tour[i + 1] - 1];

        if (strcmp(edge_weight_type, "EUC_2D") == 0) {
            total_distance += euclidean_distance(p1, p2);
        } else {
            total_distance += att_distance(p1, p2);
        }
    }

    const Point& start = points[tour[0] - 1];
    const Point& end = points[tour.back() - 1];

    if (strcmp(edge_weight_type, "EUC_2D") == 0) {
        total_distance += euclidean_distance(start, end);
    } else {
        total_distance += att_distance(start, end);
    }

    cout << "Total distance of the tour: " << total_distance << endl;

    return 0;
}
