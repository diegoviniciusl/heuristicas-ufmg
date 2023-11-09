#include <algorithm>
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

double compute_distance(const Point& a, const Point& b, const char* distance_type) {
    if (strcmp(distance_type, "EUC_2D") == 0) {
        return euclidean_distance(a, b);
    } else {
        return att_distance(a, b);
    }
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

            dist = compute_distance(points[i], current_point, distance_type) + compute_distance(points[i], next_point, distance_type) - compute_distance(current_point, next_point, distance_type);

            if (dist < min_dist) {
                min_dist = dist;
                insert_position = j;
            }
        }

        tour.insert(tour.begin() + insert_position + 1, points[i].id);
    }

    return tour;
}

bool two_opt(vector<int>& tour, const vector<Point>& points, const char* distance_type) {
    bool improved = false;

    for (size_t i = 1; i < tour.size() - 1; ++i) {
        for (size_t j = i + 1; j < tour.size(); ++j) {
            double old_distance = compute_distance(points[tour[i - 1] - 1], points[tour[i] - 1], distance_type) + compute_distance(points[tour[j] - 1], points[tour[j % tour.size()] - 1], distance_type);

            double new_distance = compute_distance(points[tour[i - 1] - 1], points[tour[j] - 1], distance_type) + compute_distance(points[tour[i] - 1], points[tour[j % tour.size()] - 1], distance_type);

            if (new_distance < old_distance) {
                reverse(tour.begin() + i, tour.begin() + j + 1);
                improved = true;
            }
        }
    }

    return improved;
}

bool three_opt(vector<int>& tour, const vector<Point>& points, const char* distance_type) {
    bool improved = false;
    for (size_t i = 0; i < tour.size() - 2; i++) {
        for (size_t j = i + 1; j < tour.size() - 1; j++) {
            for (size_t k = j + 1; k < tour.size(); k++) {
                // Compute the savings
                double delta = -compute_distance(points[tour[i] - 1], points[tour[i + 1] - 1], distance_type) - compute_distance(points[tour[j] - 1], points[tour[j + 1] - 1], distance_type) - compute_distance(points[tour[k] - 1], points[tour[(k + 1) % tour.size()] - 1], distance_type) + compute_distance(points[tour[i] - 1], points[tour[j] - 1], distance_type) + compute_distance(points[tour[i + 1] - 1], points[tour[k] - 1], distance_type) + compute_distance(points[tour[j + 1] - 1], points[tour[(k + 1) % tour.size()] - 1], distance_type);

                if (delta < 0) {
                    improved = true;
                    vector<int> new_tour;
                    // First segment
                    for (size_t l = 0; l <= i; l++) {
                        new_tour.push_back(tour[l]);
                    }
                    // Second segment reversed
                    for (size_t l = j; l > i; l--) {
                        new_tour.push_back(tour[l]);
                    }
                    // Third segment reversed
                    for (size_t l = k; l > j; l--) {
                        new_tour.push_back(tour[l]);
                    }
                    for (size_t l = k + 1; l < tour.size(); l++) {
                        new_tour.push_back(tour[l]);
                    }
                    tour = new_tour;
                }
            }
        }
    }
    return improved;
}

double calculate_total_distance(const vector<int>& tour, const vector<Point>& points, const char* distance_type) {
    double total_distance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        total_distance += compute_distance(points[tour[i] - 1], points[tour[i + 1] - 1], distance_type);
    }
    total_distance += compute_distance(points[tour.back() - 1], points[tour[0] - 1], distance_type);
    return total_distance;
}

vector<int> simulated_annealing(vector<int>& tour, const vector<Point>& points, const char* distance_type) {
    double T = 10000.0;    // Initial temperature
    double Tmin = 1.0;     // Minimum temperature
    double alpha = 0.99;   // Cooling rate
    int inner_loop = 100;  // Number of iterations at each temperature level

    srand(time(NULL));  // Seed for randomness
    double current_distance = calculate_total_distance(tour, points, distance_type);

    while (T > Tmin) {
        for (int i = 0; i < inner_loop; ++i) {
            // Create a new neighbor tour
            vector<int> new_tour = tour;
            int a = rand() % tour.size();
            int b = (a + 1 + rand() % (tour.size() - 2)) % tour.size();
            if (a > b) swap(a, b);
            reverse(new_tour.begin() + a, new_tour.begin() + b);  // This is a 2-opt move

            double new_distance = calculate_total_distance(new_tour, points, distance_type);

            // Decide if we should move to the new tour
            double ap = exp((current_distance - new_distance) / T);
            if (ap > ((double)rand() / RAND_MAX)) {
                tour = new_tour;
                current_distance = new_distance;
            }
        }

        T *= alpha;  // Cool down the temperature
    }

    return tour;  // Return the best found solution
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

    vector<int> initial_tour = nearest_insertion(points, edge_weight_type);
    vector<int> sa_tour = simulated_annealing(initial_tour, points, edge_weight_type);

    double total_distance = calculate_total_distance(sa_tour, points, edge_weight_type);
    cout << "Total distance of the tour: " << total_distance << endl;

    return 0;
}
