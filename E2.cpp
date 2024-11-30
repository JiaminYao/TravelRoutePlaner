#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <tuple>

using namespace std;

struct PathInfo {
    vector<int> path;
    double cost;
    double time;
};

// Function to load city coordinates from the CSV file
vector<cv::Point> loadCityCoordinates(const string& coor_file, vector<string>& city_names) {
    vector<cv::Point> city_list;
    ifstream infile(coor_file);
    if (!infile.is_open()) {
        cerr << "Error: Could not open the coordinates file." << endl;
        return city_list;
    }

    string line;
    getline(infile, line);  // Skip the header line if it exists

    while (getline(infile, line)) {
        stringstream ss(line);
        string name, x_str, y_str;
        int x, y;

        // Extract city name, x-coordinate, and y-coordinate
        getline(ss, name, ',');
        getline(ss, x_str, ',');
        getline(ss, y_str, ',');

        try {
            x = stoi(x_str);
            y = stoi(y_str);
            city_list.emplace_back(x, y);
            city_names.push_back(name);
        } catch (invalid_argument&) 
        {
            cerr << "Error: Invalid numeric value in coordinates for city " << name << endl;
            continue;
        }
    }

    infile.close();
    return city_list;
}

// Function to load travel costs and times
unordered_map<string, vector<double>> loadTravelData(const string& file) {
    unordered_map<string, vector<double>> travel_data;
    ifstream infile(file);
    if (!infile.is_open()) {
        cerr << "Error: Could not open the travel data file." << endl;
        return travel_data;
    }

    string line;
    getline(infile, line);  // Skip header

    while (getline(infile, line)) {
        stringstream ss(line);
        string from, to, reg_cost_str, reg_time_str, hw_cost_str, hw_time_str;
        double reg_cost, reg_time, hw_cost, hw_time;

        // Read data with commas
        getline(ss, from, ',');
        getline(ss, to, ',');
        getline(ss, reg_cost_str, ',');
        getline(ss, reg_time_str, ',');
        getline(ss, hw_cost_str, ',');
        getline(ss, hw_time_str, ',');

        // Convert strings to double
        try {
            reg_cost = stod(reg_cost_str);
            reg_time = stod(reg_time_str);
            hw_cost = stod(hw_cost_str);
            hw_time = stod(hw_time_str);

            string key = from + "_" + to;
            travel_data[key] = {reg_cost, reg_time, hw_cost, hw_time};
        } catch (const invalid_argument& e) {
            cerr << "Error: Invalid numeric data in line: " << line << endl;
            continue;
        }
    }
    infile.close();
    return travel_data;
}

// Function to validate the travel data
void validateTravelData(const vector<string>& city_names, const unordered_map<string, vector<double>>& travel_data) {
    for (int i = 0; i < city_names.size(); ++i) {
        for (int j = 0; j < city_names.size(); ++j) {
            if (i != j) {
                string key = city_names[i] + "_" + city_names[j];
                string reverse_key = city_names[j] + "_" + city_names[i];
                if (travel_data.find(key) == travel_data.end() && travel_data.find(reverse_key) == travel_data.end()) {
                    cerr << "Missing travel data for: " << city_names[i] << " <-> " << city_names[j] << endl;
                }
            }
        }
    }
}

// Function to calculate travel cost or time between two cities
double getTravelMetric(const unordered_map<string, vector<double>>& travel_data, const vector<string>& city_names, int city1, int city2, bool highway, bool cost_metric) {
    string key = city_names[city1] + "_" + city_names[city2];
    string reverse_key = city_names[city2] + "_" + city_names[city1];
    int index = highway ? (cost_metric ? 2 : 3) : (cost_metric ? 0 : 1);

    if (travel_data.find(key) != travel_data.end()) {
        return travel_data.at(key)[index];
    } else if (travel_data.find(reverse_key) != travel_data.end()) {
        return travel_data.at(reverse_key)[index];
    }
    return numeric_limits<double>::max(); // No data, return large value
}

// Function to implement the Greedy algorithm to find the min cost and min time path
std::tuple<vector<int>, vector<int>, vector<int>, vector<int>> drivingGreedyPath(
    const vector<cv::Point>& city_list,
    const unordered_map<string, vector<double>>& travel_data,
    const vector<string>& city_names,
    int start_city) {

    int n = city_list.size();
    if (n == 0) return {{}, {}, {}, {}};

    // Initialize paths and visited flags for each path type
    vector<int> regular_cost_path, regular_time_path, highway_cost_path, highway_time_path;
    vector<bool> visited_regular_cost(n, false), visited_regular_time(n, false);
    vector<bool> visited_highway_cost(n, false), visited_highway_time(n, false);

    int current_city_regular_cost = start_city, current_city_regular_time = start_city;
    int current_city_highway_cost = start_city, current_city_highway_time = start_city;

    // Start paths with the initial city and mark it as visited
    regular_cost_path.push_back(current_city_regular_cost);
    regular_time_path.push_back(current_city_regular_time);
    highway_cost_path.push_back(current_city_highway_cost);
    highway_time_path.push_back(current_city_highway_time);

    visited_regular_cost[current_city_regular_cost] = true;
    visited_regular_time[current_city_regular_time] = true;
    visited_highway_cost[current_city_highway_cost] = true;
    visited_highway_time[current_city_highway_time] = true;

    for (int i = 1; i < n; i++) {
        // Initialize minimum values and next city selections for each path type
        double min_reg_cost = numeric_limits<double>::max();
        double min_reg_time = numeric_limits<double>::max();
        double min_highway_cost = numeric_limits<double>::max();
        double min_highway_time = numeric_limits<double>::max();

        int next_city_reg_cost = -1, next_city_reg_time = -1;
        int next_city_highway_cost = -1, next_city_highway_time = -1;

        // Loop through each city to determine the next city for each path type
        for (int j = 0; j < n; j++) {
            if (!visited_regular_cost[j]) {
                double reg_cost = getTravelMetric(travel_data, city_names, current_city_regular_cost, j, false, true);
                if (reg_cost < min_reg_cost) {
                    min_reg_cost = reg_cost;
                    next_city_reg_cost = j;
                }
            }

            if (!visited_regular_time[j]) {
                double reg_time = getTravelMetric(travel_data, city_names, current_city_regular_time, j, false, false);
                if (reg_time < min_reg_time) {
                    min_reg_time = reg_time;
                    next_city_reg_time = j;
                }
            }

            if (!visited_highway_cost[j]) {
                double hw_cost = getTravelMetric(travel_data, city_names, current_city_highway_cost, j, true, true);
                if (hw_cost < min_highway_cost) {
                    min_highway_cost = hw_cost;
                    next_city_highway_cost = j;
                }
            }

            if (!visited_highway_time[j]) {
                double hw_time = getTravelMetric(travel_data, city_names, current_city_highway_time, j, true, false);
                if (hw_time < min_highway_time) {
                    min_highway_time = hw_time;
                    next_city_highway_time = j;
                }
            }
        }

        // Update paths and mark cities as visited for each path type
        if (next_city_reg_cost != -1) {
            regular_cost_path.push_back(next_city_reg_cost);
            visited_regular_cost[next_city_reg_cost] = true;
            current_city_regular_cost = next_city_reg_cost;
        }

        if (next_city_reg_time != -1) {
            regular_time_path.push_back(next_city_reg_time);
            visited_regular_time[next_city_reg_time] = true;
            current_city_regular_time = next_city_reg_time;
        }

        if (next_city_highway_cost != -1) {
            highway_cost_path.push_back(next_city_highway_cost);
            visited_highway_cost[next_city_highway_cost] = true;
            current_city_highway_cost = next_city_highway_cost;
        }

        if (next_city_highway_time != -1) {
            highway_time_path.push_back(next_city_highway_time);
            visited_highway_time[next_city_highway_time] = true;
            current_city_highway_time = next_city_highway_time;
        }
    }

    return {regular_cost_path, regular_time_path, highway_cost_path, highway_time_path};
}

// Function to draw the Shortest Path using Greedy Algorithm
void drawDrivingGreedyPath(cv::Mat& base_img, const vector<cv::Point>& city_coords, const vector<int>& path, const string& output_file) {
    cv::Mat img = base_img.clone();
    for (size_t i = 0; i < path.size(); ++i) {
        int greedy_index = path[i];
        cv::circle(img, city_coords[greedy_index], 10, cv::Scalar(0, 0, 255), -1);
        cv::putText(img, to_string(i), city_coords[greedy_index], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    for (size_t i = 0; i < path.size() - 1; i++) {
        cv::line(img, city_coords[path[i]], city_coords[path[i + 1]], cv::Scalar(30, 150, 80), 4);
    }
    cv::line(img, city_coords[path.back()], city_coords[path[0]], cv::Scalar(30, 150, 80), 4);
    cv::imwrite(output_file, img);
}



// Recursive function to find the optimal path by Divide-and-Conquer without skipping cities
PathInfo divideAndConquerPath(int start, int end, const unordered_map<string, vector<double>>& travel_data, const vector<string>& city_names, bool highway, bool cost_metric, unordered_set<int>& visited) {
    // Base case
    if (start > end) return {{}, 0, 0};
    if (start == end) {
        visited.insert(start);
        return {{start}, 0, 0};
    }

    // Collect all cities in the range and their metrics
    vector<pair<int, double>> cities;
    for (int i = start; i <= end; ++i) {
        double metric_value = getTravelMetric(travel_data, city_names, start, i, highway, cost_metric);
        cities.emplace_back(i, metric_value);
    }

    // Sort cities by the selected metric (cost or time)
    sort(cities.begin(), cities.end(), [cost_metric](const pair<int, double>& a, const pair<int, double>& b) {
        return a.second < b.second;
    });

    // Extract sorted city indices
    vector<int> sorted_indices;
    for (const auto& city : cities) {
        sorted_indices.push_back(city.first);
    }

    // Divide
    int mid = sorted_indices.size() / 2;
    PathInfo left_path = divideAndConquerPath(sorted_indices.front(), sorted_indices[mid - 1], travel_data, city_names, highway, cost_metric, visited);
    PathInfo right_path = divideAndConquerPath(sorted_indices[mid], sorted_indices.back(), travel_data, city_names, highway, cost_metric, visited);

    // Calculate the cost/time to connect the last city in the left path to the first city in the right path
    double min_bridge_cost = numeric_limits<double>::max();
    double min_bridge_time = numeric_limits<double>::max();
    int best_left = -1, best_right = -1;

    for (int left_city : left_path.path) {
        for (int right_city : right_path.path) {
            if (visited.find(right_city) != visited.end()) continue;

            double bridge_cost = getTravelMetric(travel_data, city_names, left_city, right_city, highway, true);
            double bridge_time = getTravelMetric(travel_data, city_names, left_city, right_city, highway, false);

            if ((cost_metric && bridge_cost < min_bridge_cost) || (!cost_metric && bridge_time < min_bridge_time)) {
                min_bridge_cost = bridge_cost;
                min_bridge_time = bridge_time;
                best_left = left_city;
                best_right = right_city;
            }
        }
    }

    // Merge paths with the optimal bridge
    PathInfo result;
    result.path = left_path.path;
    if (best_right != -1 && visited.find(best_right) == visited.end()) {
        result.path.push_back(best_right);  // Add the bridge city
        visited.insert(best_right);        // Mark as visited
    }
    result.path.insert(result.path.end(), right_path.path.begin(), right_path.path.end());

    result.cost = left_path.cost + right_path.cost + min_bridge_cost;
    result.time = left_path.time + right_path.time + min_bridge_time;

    // Ensure all cities in the range are visited
    for (int i = start; i <= end; ++i) {
        if (visited.find(i) == visited.end()) {
            result.path.push_back(i);
            visited.insert(i);
        }
    }

    return result;
}

// Wrapper function to generate paths for all metrics using Divide-and-Conquer
std::tuple<vector<int>, vector<int>, vector<int>, vector<int>> drivingDCPath(
    const vector<cv::Point>& city_list,
    const unordered_map<string, vector<double>>& travel_data,
    const vector<string>& city_names,
    int start_city) {

    int n = city_list.size();
    if (n == 0) return {{}, {}, {}, {}};

    unordered_set<int> visited;

    PathInfo reg_cost_path = divideAndConquerPath(0, n - 1, travel_data, city_names, false, true, visited);
    visited.clear();
    PathInfo reg_time_path = divideAndConquerPath(0, n - 1, travel_data, city_names, false, false, visited);
    visited.clear();
    PathInfo hw_cost_path = divideAndConquerPath(0, n - 1, travel_data, city_names, true, true, visited);
    visited.clear();
    PathInfo hw_time_path = divideAndConquerPath(0, n - 1, travel_data, city_names, true, false, visited);

    // Ensure the start city is at the beginning of each path
    auto adjustPathStart = [&](vector<int>& path) {
        auto it = find(path.begin(), path.end(), start_city);
        if (it != path.end()) {
            rotate(path.begin(), it, path.end());
        }
    };

    adjustPathStart(reg_cost_path.path);
    adjustPathStart(reg_time_path.path);
    adjustPathStart(hw_cost_path.path);
    adjustPathStart(hw_time_path.path);

    return {reg_cost_path.path, reg_time_path.path, hw_cost_path.path, hw_time_path.path};
}

void drawDrivingDCPath(cv::Mat& base_img, const vector<cv::Point>& city_coords, const vector<int>& path, const string& output_file) {
    cv::Mat img = base_img.clone();
    for (size_t i = 0; i < path.size(); ++i) {
        int dc_index = path[i];
        cv::circle(img, city_coords[dc_index], 10, cv::Scalar(0, 0, 255), -1);
        cv::putText(img, to_string(i), city_coords[dc_index], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    for (size_t i = 0; i < path.size() - 1; ++i) {
        cv::line(img, city_coords[path[i]], city_coords[path[i + 1]], cv::Scalar(170, 30, 220), 4);
    }
    if (!path.empty()) {
        cv::line(img, city_coords[path.back()], city_coords[path[0]], cv::Scalar(170, 30, 220), 4);
    }
    cv::imwrite(output_file, img);
}

// Precompute travel metrics for faster DP calculation
vector<vector<double>> precomputeTravelMetrics(
    const unordered_map<string, vector<double>>& travel_data,
    const vector<string>& city_names, bool highway, bool cost_metric) {

    int n = city_names.size();
    vector<vector<double>> travel_metric(n, vector<double>(n, numeric_limits<double>::max()));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                string key = city_names[i] + "_" + city_names[j];
                string reverse_key = city_names[j] + "_" + city_names[i];
                travel_metric[i][j] = travel_data.count(key) ? travel_data.at(key)[cost_metric ? 0 : 1] :
                                        (travel_data.count(reverse_key) ? travel_data.at(reverse_key)[cost_metric ? 0 : 1] : numeric_limits<double>::max());
            }
        }
    }

    return travel_metric;
}

// DP function to solve TSP
PathInfo dynamicProgrammingTSP(const vector<vector<double>>& travel_metric) {
    int n = travel_metric.size();
    if (n == 0) return {{}, 0.0, 0.0};

    // DP table to store minimum cost of visiting cities with a specific bitmask
    vector<vector<double>> dp_cost(1 << n, vector<double>(n, -1.0));
    vector<vector<int>> parent_cost(1 << n, vector<int>(n, -1));

    // Recursive helper function
    function<double(int, int)> dpHelper = [&](int current_city, int visited_mask) {
        if (visited_mask == (1 << n) - 1) {
            // Base case: all cities visited, return to the starting city
            return travel_metric[current_city][0];
        }

        if (dp_cost[visited_mask][current_city] != -1.0) {
            // Return memoized result if available
            return dp_cost[visited_mask][current_city];
        }

        double min_cost = numeric_limits<double>::max();
        int next_city = -1;

        // Explore all unvisited cities
        for (int next = 0; next < n; ++next) {
            if (!(visited_mask & (1 << next))) { // If the city is not visited
                double cost = travel_metric[current_city][next] + dpHelper(next, visited_mask | (1 << next));
                if (cost < min_cost) {
                    min_cost = cost;
                    next_city = next;
                }
            }
        }

        dp_cost[visited_mask][current_city] = min_cost;
        parent_cost[visited_mask][current_city] = next_city;
        return min_cost;
    };

    // Start the DP recursion from city 0 with only it visited
    dpHelper(0, 1);

    // Reconstruct the path
    vector<int> path;
    int mask = 1, current_city = 0;

    while (current_city != -1) {
        path.push_back(current_city);
        int next_city = parent_cost[mask][current_city];
        mask |= (1 << next_city);
        current_city = next_city;
    }

    // Add the starting city at the end to complete the cycle
    path.push_back(0);

    // Compute total cost of the path
    double total_cost = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        total_cost += travel_metric[path[i]][path[i + 1]];
    }

    return {path, total_cost, 0.0};
}

PathInfo dynamicProgrammingTSPWithStart(const vector<vector<double>>& travel_metric, int start_city) {
    int n = travel_metric.size();
    if (n == 0) return {{}, 0.0, 0.0};

    // DP table to store the minimum cost of visiting cities with a specific bitmask
    vector<vector<double>> dp_cost(1 << n, vector<double>(n, -1.0));
    vector<vector<int>> parent_cost(1 << n, vector<int>(n, -1));

    // Recursive helper function
    function<double(int, int)> dpHelper = [&](int current_city, int visited_mask) {
        if (visited_mask == (1 << n) - 1) {
            // Base case: all cities visited, return to the starting city
            return travel_metric[current_city][start_city];
        }

        if (dp_cost[visited_mask][current_city] != -1.0) {
            // Return memoized result if available
            return dp_cost[visited_mask][current_city];
        }

        double min_cost = numeric_limits<double>::max();
        int next_city = -1;

        // Explore all unvisited cities
        for (int next = 0; next < n; ++next) {
            if (!(visited_mask & (1 << next))) { // If the city is not visited
                double cost = travel_metric[current_city][next] + dpHelper(next, visited_mask | (1 << next));
                if (cost < min_cost) {
                    min_cost = cost;
                    next_city = next;
                }
            }
        }

        dp_cost[visited_mask][current_city] = min_cost;
        parent_cost[visited_mask][current_city] = next_city;
        return min_cost;
    };

    // Start the DP recursion from the specified starting city with only it visited
    dpHelper(start_city, 1 << start_city);

    // Reconstruct the path
    vector<int> path;
    int mask = 1 << start_city, current_city = start_city;

    while (current_city != -1) {
        path.push_back(current_city);
        int next_city = parent_cost[mask][current_city];
        mask |= (1 << next_city);
        current_city = next_city;
    }

    // Add the starting city at the end to complete the cycle
    path.push_back(start_city);

    // Compute total cost of the path
    double total_cost = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        total_cost += travel_metric[path[i]][path[i + 1]];
    }

    return {path, total_cost, 0.0};
}


// Wrapper function for DP paths
std::tuple<vector<int>, vector<int>, vector<int>, vector<int>> drivingDPPath(
    const vector<cv::Point>& city_list,
    const unordered_map<string, vector<double>>& travel_data,
    const vector<string>& city_names,
    int start_city) {

    int n = city_list.size();
    if (n == 0) return {{}, {}, {}, {}};

    // Precompute travel metrics
    vector<vector<double>> reg_cost_metric = precomputeTravelMetrics(travel_data, city_names, false, true);
    vector<vector<double>> reg_time_metric = precomputeTravelMetrics(travel_data, city_names, false, false);
    vector<vector<double>> hw_cost_metric = precomputeTravelMetrics(travel_data, city_names, true, true);
    vector<vector<double>> hw_time_metric = precomputeTravelMetrics(travel_data, city_names, true, false);

    // Calculate DP paths starting from the selected city
    PathInfo reg_cost_path = dynamicProgrammingTSPWithStart(reg_cost_metric, start_city);
    PathInfo reg_time_path = dynamicProgrammingTSPWithStart(reg_time_metric, start_city);
    PathInfo hw_cost_path = dynamicProgrammingTSPWithStart(hw_cost_metric, start_city);
    PathInfo hw_time_path = dynamicProgrammingTSPWithStart(hw_time_metric, start_city);

    return {reg_cost_path.path, reg_time_path.path, hw_cost_path.path, hw_time_path.path};
}


void drawDrivingDPPath(cv::Mat& base_img, const vector<cv::Point>& city_coords, const vector<int>& path, const string& output_file) {
    cv::Mat img = base_img.clone();

    // Draw cities and indices along the path
    for (size_t i = 0; i < path.size() -1; ++i) {
        int index = path[i];
        cv::circle(img, city_coords[index], 10, cv::Scalar(0, 0, 255), -1); // Green for cities
        cv::putText(img, to_string(i), city_coords[index], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw connections between cities
    for (size_t i = 0; i < path.size() - 1; ++i) {
        cv::line(img, city_coords[path[i]], city_coords[path[i + 1]], cv::Scalar(0, 255, 255), 4); // Red for paths
    }

    // Connect the last city back to the first
    if (!path.empty()) {
        cv::line(img, city_coords[path.back()], city_coords[path[0]], cv::Scalar(0, 255, 255), 4); // Closing the path
    }

    // Save the visualization
    cv::imwrite(output_file, img);
}

// Function to calculate the total cost or time for a given path
double calculateTotalMetric(const vector<int>& path, const unordered_map<string, vector<double>>& travel_data, const vector<string>& city_names, bool cost_metric = true, bool highway = false) {
    double total_metric = 0.0;

    for (size_t i = 0; i < path.size() - 1; ++i) {
        if (path[i] == path[i + 1]) { // Skip self-loops
            continue;
        }

        string key = city_names[path[i]] + "_" + city_names[path[i + 1]];
        string reverse_key = city_names[path[i + 1]] + "_" + city_names[path[i]];

        int index = highway ? (cost_metric ? 2 : 3) : (cost_metric ? 0 : 1); // Choose the correct index based on highway and metric type

        if (travel_data.find(key) != travel_data.end()) {
            total_metric += travel_data.at(key)[index];
        } else if (travel_data.find(reverse_key) != travel_data.end()) {
            total_metric += travel_data.at(reverse_key)[index];
        } else {
            cerr << "Warning: No travel data found for route " << key << " or " << reverse_key << endl;
        }
    }

    return total_metric;
}

void writeCityOrderToCSV(const string& filename, const vector<int>& path, const vector<string>& city_names, const vector<cv::Point>& city_coords, int total_cost, int total_time) {
    ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << "Index,City Name,Coordinate (x,y)" << endl;

        for (size_t new_index = 0; new_index < path.size(); ++new_index) {
            // Skip the last city if it's a duplicate of the first city
            if (new_index == path.size() - 1 && path[new_index] == path[0]) {
                continue;
            }

            int original_index = path[new_index];
            outfile << new_index << "," << city_names[original_index]
                    << ",(" << city_coords[original_index].x 
                    << "," << city_coords[original_index].y << ")" << endl;
        }
        outfile.close();
    } else {
        cerr << "Error: Could not open file " << filename << " for writing." << endl;
    }
}

int main() {
    vector<string> city_names;
    vector<cv::Point> city_coords = loadCityCoordinates("./Dataset/Dataset_Coordinate.csv", city_names);
    if (city_coords.empty()) return -1;

    auto travel_data = loadTravelData("./Dataset/Dataset_Driving.csv");
    validateTravelData(city_names, travel_data);
    cv::Mat base_img = cv::imread("./Image/Europe.png");
    if (base_img.empty()) return -1;

    cout << "Available cities:" << endl;
    for (size_t i = 0; i < city_names.size(); ++i) {
        cout << i << ": " << city_names[i] << " (" << city_coords[i].x << ", " << city_coords[i].y << ")" << endl;
    }

    int start_city_index;
    cout << "Enter the index of the starting city: ";
    cin >> start_city_index;

    // Greedy paths
    auto [reg_cost_path_greedy, reg_time_path_greedy, hw_cost_path_greedy, hw_time_path_greedy] =
        drivingGreedyPath(city_coords, travel_data, city_names, start_city_index);

    // Divide-and-Conquer paths
    auto [reg_cost_path_dc, reg_time_path_dc, hw_cost_path_dc, hw_time_path_dc] =
        drivingDCPath(city_coords, travel_data, city_names, start_city_index);

    // Dynamic Programming paths
    auto [reg_cost_path_dp, reg_time_path_dp, hw_cost_path_dp, hw_time_path_dp] =
        drivingDPPath(city_coords, travel_data, city_names, start_city_index);

    // Define path information for each type
    struct PathInfo {
        vector<int> path;
        bool highway;
        string path_name;
        string image_filename;
        string csv_filename;
    };

    vector<PathInfo> greedy_paths = {
        {reg_cost_path_greedy, false, "Regular Road Min Cost Path", "./Image/E2Driving_Greedy_MinRegularCost_Path.png", "./Table/E2Driving_Greedy_MinRegularCost_Order.csv"},
        {reg_time_path_greedy, false, "Regular Road Min Time Path", "./Image/E2Driving_Greedy_MinRegularTime_Path.png", "./Table/E2Driving_Greedy_MinRegularTime_Order.csv"},
        {hw_cost_path_greedy, true, "Highway Min Cost Path", "./Image/E2Driving_Greedy_MinHighwayCost_Path.png", "./Table/E2Driving_Greedy_MinHighwayCost_Order.csv"},
        {hw_time_path_greedy, true, "Highway Min Time Path", "./Image/E2Driving_Greedy_MinHighwayTime_Path.png", "./Table/E2Driving_Greedy_MinHighwayTime_Order.csv"}
    };

    vector<PathInfo> dc_paths = {
        {reg_cost_path_dc, false, "Regular Road Min Cost Path", "./Image/E2Driving_DC_MinRegularCost_Path.png", "./Table/E2Driving_DC_MinRegularCost_Order.csv"},
        {reg_time_path_dc, false, "Regular Road Min Time Path", "./Image/E2Driving_DC_MinRegularTime_Path.png", "./Table/E2Driving_DC_MinRegularTime_Order.csv"},
        {hw_cost_path_dc, true, "Highway Min Cost Path", "./Image/E2Driving_DC_MinHighwayCost_Path.png", "./Table/E2Driving_DC_MinHighwayCost_Order.csv"},
        {hw_time_path_dc, true, "Highway Min Time Path", "./Image/E2Driving_DC_MinHighwayTime_Path.png", "./Table/E2Driving_DC_MinHighwayTime_Order.csv"}
    };

    vector<PathInfo> dp_paths = {
        {reg_cost_path_dp, false, "Regular Road Min Cost Path", "./Image/E2Driving_DP_MinRegularCost_Path.png", "./Table/E2Driving_DP_MinRegularCost_Order.csv"},
        {reg_time_path_dp, false, "Regular Road Min Time Path", "./Image/E2Driving_DP_MinRegularTime_Path.png", "./Table/E2Driving_DP_MinRegularTime_Order.csv"},
        {hw_cost_path_dp, true, "Highway Min Cost Path", "./Image/E2Driving_DP_MinHighwayCost_Path.png", "./Table/E2Driving_DP_MinHighwayCost_Order.csv"},
        {hw_time_path_dp, true, "Highway Min Time Path", "./Image/E2Driving_DP_MinHighwayTime_Path.png", "./Table/E2Driving_DP_MinHighwayTime_Order.csv"}
    };

    // Function to process paths for any algorithm
    auto processPaths = [&](const vector<PathInfo>& paths, const string& algorithm) {
        for (const auto& path_info : paths) {
            double total_cost = calculateTotalMetric(path_info.path, travel_data, city_names, true, path_info.highway);
            double total_time = calculateTotalMetric(path_info.path, travel_data, city_names, false, path_info.highway);

            cout << algorithm << " " << path_info.path_name
                    << " Total Cost: $" << total_cost << ", Total Time: " << total_time << " mins" << endl;

            cv::Mat img = cv::imread("./Image/Europe.png");
            if (img.empty()) {
                cerr << "Error: Could not load image for path visualization." << endl;
                continue;
            }

            if (algorithm == "Greedy") {
                drawDrivingGreedyPath(img, city_coords, path_info.path, path_info.image_filename);
            } else if (algorithm == "DC") {
                drawDrivingDCPath(img, city_coords, path_info.path, path_info.image_filename);
            } else if (algorithm == "DP") {
                drawDrivingDPPath(img, city_coords, path_info.path, path_info.image_filename);
            }

            writeCityOrderToCSV(path_info.csv_filename, path_info.path, city_names, city_coords, total_cost, total_time);
        }
    };

    // Process and save paths for all algorithms
    processPaths(greedy_paths, "Greedy");
    processPaths(dc_paths, "DC");
    processPaths(dp_paths, "DP");

    // Open the saved images using the system's default image viewer
    std::vector<PathInfo> all_paths = greedy_paths;
    all_paths.insert(all_paths.end(), dc_paths.begin(), dc_paths.end()); // Add Divide & Conquer paths
    all_paths.insert(all_paths.end(), dp_paths.begin(), dp_paths.end()); // Add Dynamic Programming paths

    for (const auto& path_info : all_paths) {
#ifdef _WIN32
        system(("start " + path_info.image_filename).c_str());
#elif __APPLE__
        system(("open " + path_info.image_filename).c_str());
#elif __linux__
        system(("xdg-open " + path_info.image_filename).c_str());
#endif
    }

    return 0;
}

