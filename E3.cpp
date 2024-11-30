#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <limits>
#include <numeric>
#include <functional>

using namespace std;

// Divide and Conquer Algorithm to find the optimal path
struct PathInfo {
    vector<int> path;
    double total_cost;
    double total_time;
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

// Function to load flight data from CSV
unordered_map<string, vector<double>> loadFlightData(const string& file, vector<string>& city_names) {
    unordered_map<string, vector<double>> flight_data;
    ifstream infile(file);
    if (!infile.is_open()) {
        cerr << "Error: Could not open the flight data file." << endl;
        return flight_data;
    }

    string line;
    getline(infile, line);  // Skip header

    while (getline(infile, line)) {
        stringstream ss(line);
        string from, to, cost_str, time_str;
        double cost, time;

        // Parse the line
        getline(ss, from, ',');
        getline(ss, to, ',');
        getline(ss, cost_str, ',');
        getline(ss, time_str, ',');

        try {
            cost = stod(cost_str);
            time = stod(time_str);

            string key = from + "_" + to;
            flight_data[key] = {cost, time};

            // Add unique city names
            if (find(city_names.begin(), city_names.end(), from) == city_names.end()) {
                city_names.push_back(from);
            }
            if (find(city_names.begin(), city_names.end(), to) == city_names.end()) {
                city_names.push_back(to);
            }

        } catch (const invalid_argument& e) {
            cerr << "Error: Invalid numeric data in line: " << line << endl;
            continue;
        }
    }
    infile.close();
    return flight_data;
}

// Greedy Algorithm to find minimum-cost and minimum-time paths
tuple<vector<int>, vector<int>> flightGreedyPath(
    const vector<cv::Point>& city_coords,
    const unordered_map<string, vector<double>>& flight_data,
    const vector<string>& city_names,
    int start_city) {

    int n = city_names.size();
    if (n == 0) return {{}, {}};

    vector<int> min_cost_path, min_time_path;
    vector<bool> visited_cost(n, false), visited_time(n, false);

    int current_city_cost = start_city, current_city_time = start_city;

    // Start paths with the selected city
    min_cost_path.push_back(current_city_cost);
    min_time_path.push_back(current_city_time);

    visited_cost[current_city_cost] = true;
    visited_time[current_city_time] = true;

    for (int i = 1; i < n; ++i) {
        double min_cost = numeric_limits<double>::max();
        double min_time = numeric_limits<double>::max();

        int next_city_cost = -1, next_city_time = -1;

        for (int j = 0; j < n; ++j) {
            if (!visited_cost[j]) {
                string key = city_names[current_city_cost] + "_" + city_names[j];
                string reverse_key = city_names[j] + "_" + city_names[current_city_cost];

                double cost = flight_data.count(key) ? flight_data.at(key)[0] :
                                (flight_data.count(reverse_key) ? flight_data.at(reverse_key)[0] : numeric_limits<double>::max());

                if (cost < min_cost) {
                    min_cost = cost;
                    next_city_cost = j;
                }
            }

            if (!visited_time[j]) {
                string key = city_names[current_city_time] + "_" + city_names[j];
                string reverse_key = city_names[j] + "_" + city_names[current_city_time];

                double time = flight_data.count(key) ? flight_data.at(key)[1] :
                                (flight_data.count(reverse_key) ? flight_data.at(reverse_key)[1] : numeric_limits<double>::max());

                if (time < min_time) {
                    min_time = time;
                    next_city_time = j;
                }
            }
        }

        if (next_city_cost != -1) {
            min_cost_path.push_back(next_city_cost);
            visited_cost[next_city_cost] = true;
            current_city_cost = next_city_cost;
        }

        if (next_city_time != -1) {
            min_time_path.push_back(next_city_time);
            visited_time[next_city_time] = true;
            current_city_time = next_city_time;
        }
    }

    return {min_cost_path, min_time_path};
}


// Function to draw the flight paths using Greedy Algorithm
void drawFlightGreedyPath(cv::Mat& base_img, const vector<cv::Point>& city_coords, const vector<int>& path, const string& output_file) {
    cv::Mat img = base_img.clone();

    // Draw cities and indices along the path
    for (size_t i = 0; i < path.size(); ++i) {
        int index = path[i];
        cv::circle(img, city_coords[index], 10, cv::Scalar(0, 0, 255), -1); // Red for cities
        cv::putText(img, to_string(i), city_coords[index], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw connections between cities
    for (size_t i = 0; i < path.size() - 1; ++i) {
        cv::line(img, city_coords[path[i]], city_coords[path[i + 1]], cv::Scalar(30, 150, 80), 4); // Green for paths
    }

    // Connect the last city back to the first
    if (!path.empty()) {
        cv::line(img, city_coords[path.back()], city_coords[path[0]], cv::Scalar(30, 150, 80), 4);
    }

    // Save the visualization
    cv::imwrite(output_file, img);
}

// Helper function for Divide and Conquer
PathInfo divideAndConquerHelper(
    int start, int end,
    const unordered_map<string, vector<double>>& flight_data,
    const vector<string>& city_names,
    bool cost_metric) {

    // Base case: only one city
    if (start == end) {
        return {{start}, 0.0, 0.0};
    }

    // Collect city pairs and their metrics
    vector<pair<int, double>> metrics;
    for (int i = start; i <= end; ++i) {
        double metric_value = 0.0;

        if (i > start) {
            string key = city_names[start] + "_" + city_names[i];
            string reverse_key = city_names[i] + "_" + city_names[start];

            if (flight_data.count(key)) {
                metric_value = flight_data.at(key)[cost_metric ? 0 : 1];
            } else if (flight_data.count(reverse_key)) {
                metric_value = flight_data.at(reverse_key)[cost_metric ? 0 : 1];
            } else {
                metric_value = numeric_limits<double>::max();
            }
        }

        metrics.push_back({i, metric_value});
    }

    // Sort by cost or time
    sort(metrics.begin(), metrics.end(), [](const pair<int, double>& a, const pair<int, double>& b) {
        return a.second < b.second;
    });

    // Extract sorted indices
    vector<int> sorted_indices;
    for (const auto& metric : metrics) {
        sorted_indices.push_back(metric.first);
    }

    // Divide indices into left and right halves
    int mid = sorted_indices.size() / 2;
    vector<int> left_indices(sorted_indices.begin(), sorted_indices.begin() + mid);
    vector<int> right_indices(sorted_indices.begin() + mid, sorted_indices.end());

    // Recursive calls for left and right halves
    PathInfo left = divideAndConquerHelper(left_indices.front(), left_indices.back(), flight_data, city_names, cost_metric);
    PathInfo right = divideAndConquerHelper(right_indices.front(), right_indices.back(), flight_data, city_names, cost_metric);

    // Merge left and right paths
    double min_bridge_cost = numeric_limits<double>::max();
    double min_bridge_time = numeric_limits<double>::max();
    int best_left = -1, best_right = -1;

    for (int l : left.path) {
        for (int r : right.path) {
            string key = city_names[l] + "_" + city_names[r];
            string reverse_key = city_names[r] + "_" + city_names[l];

            double bridge_cost = flight_data.count(key) ? flight_data.at(key)[0] :
                                    (flight_data.count(reverse_key) ? flight_data.at(reverse_key)[0] : numeric_limits<double>::max());
            double bridge_time = flight_data.count(key) ? flight_data.at(key)[1] :
                                    (flight_data.count(reverse_key) ? flight_data.at(reverse_key)[1] : numeric_limits<double>::max());

            if (cost_metric && bridge_cost < min_bridge_cost) {
                min_bridge_cost = bridge_cost;
                min_bridge_time = bridge_time;
                best_left = l;
                best_right = r;
            } else if (!cost_metric && bridge_time < min_bridge_time) {
                min_bridge_cost = bridge_cost;
                min_bridge_time = bridge_time;
                best_left = l;
                best_right = r;
            }
        }
    }

    PathInfo result;
    result.path = left.path;
    result.path.insert(result.path.end(), right.path.begin(), right.path.end());
    result.total_cost = left.total_cost + right.total_cost + min_bridge_cost;
    result.total_time = left.total_time + right.total_time + min_bridge_time;

    return result;
}


// Wrapper function for Divide and Conquer
tuple<vector<int>, vector<int>> flightDCPath(
    const vector<cv::Point>& city_coords,
    const unordered_map<string, vector<double>>& flight_data,
    const vector<string>& city_names,
    int start_city) {

    int n = city_names.size();
    if (n == 0) return {{}, {}};

    // Adjust indices for the selected starting city
    vector<int> indices(n);
    iota(indices.begin(), indices.end(), 0);
    rotate(indices.begin(), indices.begin() + start_city, indices.end());

    // Recursive Divide and Conquer helper
    std::function<PathInfo(int, int, bool)> divideAndConquerHelper =
        [&](int start, int end, bool cost_metric) -> PathInfo {
        if (start == end) {
            return {{indices[start]}, 0.0, 0.0};
        }

        int mid = (start + end) / 2;
        PathInfo left = divideAndConquerHelper(start, mid, cost_metric);
        PathInfo right = divideAndConquerHelper(mid + 1, end, cost_metric);

        // Merge left and right paths
        double min_bridge_cost = numeric_limits<double>::max();
        double min_bridge_time = numeric_limits<double>::max();
        int best_left = -1, best_right = -1;

        for (int l : left.path) {
            for (int r : right.path) {
                string key = city_names[l] + "_" + city_names[r];
                string reverse_key = city_names[r] + "_" + city_names[l];

                double bridge_cost = flight_data.count(key) ? flight_data.at(key)[0] :
                                        (flight_data.count(reverse_key) ? flight_data.at(reverse_key)[0] : numeric_limits<double>::max());
                double bridge_time = flight_data.count(key) ? flight_data.at(key)[1] :
                                        (flight_data.count(reverse_key) ? flight_data.at(reverse_key)[1] : numeric_limits<double>::max());

                if (cost_metric && bridge_cost < min_bridge_cost) {
                    min_bridge_cost = bridge_cost;
                    min_bridge_time = bridge_time;
                    best_left = l;
                    best_right = r;
                } else if (!cost_metric && bridge_time < min_bridge_time) {
                    min_bridge_cost = bridge_cost;
                    min_bridge_time = bridge_time;
                    best_left = l;
                    best_right = r;
                }
            }
        }

        PathInfo result;
        result.path = left.path;
        result.path.insert(result.path.end(), right.path.begin(), right.path.end());
        result.total_cost = left.total_cost + right.total_cost + min_bridge_cost;
        result.total_time = left.total_time + right.total_time + min_bridge_time;

        return result;
    };

    // Find paths
    PathInfo min_cost_path_info = divideAndConquerHelper(0, n - 1, true);
    PathInfo min_time_path_info = divideAndConquerHelper(0, n - 1, false);

    return {min_cost_path_info.path, min_time_path_info.path};
}

// Function to draw the flight paths using Divide and Conquer Algorithm
void drawFlightDCPath(cv::Mat& base_img, const vector<cv::Point>& city_coords, const vector<int>& path, const string& output_file) {
    cv::Mat img = base_img.clone();

    // Draw cities and indices along the path
    for (size_t i = 0; i < path.size(); ++i) {
        int index = path[i];
        cv::circle(img, city_coords[index], 10, cv::Scalar(0, 0, 255), -1); // Blue for cities
        cv::putText(img, to_string(i), city_coords[index], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw connections between cities
    for (size_t i = 0; i < path.size() - 1; ++i) {
        cv::line(img, city_coords[path[i]], city_coords[path[i + 1]], cv::Scalar(170, 30, 220), 4); // Green for paths
    }

    // Connect the last city back to the first
    if (!path.empty()) {
        cv::line(img, city_coords[path.back()], city_coords[path[0]], cv::Scalar(170, 30, 220), 4); // Green for closing the path
    }

    // Save the visualization
    cv::imwrite(output_file, img);
}

// Dynamic Programming to find the optimal flight path
PathInfo dynamicProgrammingHelper(
    int current_city,
    int visited_mask,
    int n,
    const vector<vector<double>>& cost_matrix,
    vector<vector<double>>& dp_cost,
    vector<vector<int>>& parent_cost) {

    if (visited_mask == (1 << n) - 1) {
        // Base case: all cities visited, return to the starting city
        return {{}, cost_matrix[current_city][0], cost_matrix[current_city][0]};
    }

    if (dp_cost[visited_mask][current_city] != -1.0) {
        // Use memoized result if available
        return {{}, dp_cost[visited_mask][current_city], 0.0};
    }

    double min_cost = numeric_limits<double>::max();
    int next_city = -1;

    // Recurse for unvisited cities
    for (int next = 0; next < n; ++next) {
        if (!(visited_mask & (1 << next))) {  // If the city is not visited
            double cost = cost_matrix[current_city][next] + 
                            dynamicProgrammingHelper(next, visited_mask | (1 << next), n, cost_matrix, dp_cost, parent_cost).total_cost;

            if (cost < min_cost) {
                min_cost = cost;
                next_city = next;
            }
        }
    }

    dp_cost[visited_mask][current_city] = min_cost;
    parent_cost[visited_mask][current_city] = next_city;
    return {{}, min_cost, 0.0};
}

// Wrapper function for Dynamic Programming
tuple<vector<int>, vector<int>> flightDPPath(
    const vector<cv::Point>& city_coords,
    const unordered_map<string, vector<double>>& flight_data,
    const vector<string>& city_names,
    int start_city) {

    int n = city_names.size();
    if (n == 0) return {{}, {}};

    // Prepare cost matrices for both metrics (cost and time)
    vector<vector<double>> cost_matrix(n, vector<double>(n, numeric_limits<double>::max()));
    vector<vector<double>> time_matrix(n, vector<double>(n, numeric_limits<double>::max()));

    // Populate cost and time matrices
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                string key = city_names[i] + "_" + city_names[j];
                string reverse_key = city_names[j] + "_" + city_names[i];

                cost_matrix[i][j] = flight_data.count(key) ? flight_data.at(key)[0] :
                                    (flight_data.count(reverse_key) ? flight_data.at(reverse_key)[0] : numeric_limits<double>::max());
                time_matrix[i][j] = flight_data.count(key) ? flight_data.at(key)[1] :
                                    (flight_data.count(reverse_key) ? flight_data.at(reverse_key)[1] : numeric_limits<double>::max());
            }
        }
    }

    // DP arrays and parent tracking
    vector<vector<double>> dp_cost(1 << n, vector<double>(n, -1.0));
    vector<vector<int>> parent_cost(1 << n, vector<int>(n, -1));

    vector<vector<double>> dp_time(1 << n, vector<double>(n, -1.0));
    vector<vector<int>> parent_time(1 << n, vector<int>(n, -1));

    // Recursive helper for DP
    std::function<double(int, int, const vector<vector<double>>&, vector<vector<double>>&, vector<vector<int>>&)> dpHelper =
        [&](int current_city, int visited_mask, const vector<vector<double>>& metric_matrix,
            vector<vector<double>>& dp, vector<vector<int>>& parent) -> double {
        if (visited_mask == (1 << n) - 1) {
            return metric_matrix[current_city][start_city];
        }

        if (dp[visited_mask][current_city] != -1.0) {
            return dp[visited_mask][current_city];
        }

        double min_metric = numeric_limits<double>::max();
        int next_city = -1;

        for (int next = 0; next < n; ++next) {
            if (!(visited_mask & (1 << next))) {  // If the city is not visited
                double cost = metric_matrix[current_city][next] +
                                dpHelper(next, visited_mask | (1 << next), metric_matrix, dp, parent);

                if (cost < min_metric) {
                    min_metric = cost;
                    next_city = next;
                }
            }
        }

        dp[visited_mask][current_city] = min_metric;
        parent[visited_mask][current_city] = next_city;
        return min_metric;
    };

    // Compute minimum-cost path
    dpHelper(start_city, 1 << start_city, cost_matrix, dp_cost, parent_cost);

    // Compute minimum-time path
    dpHelper(start_city, 1 << start_city, time_matrix, dp_time, parent_time);

    // Reconstruct paths
    auto reconstructPath = [&](const vector<vector<int>>& parent) {
        vector<int> path;
        int mask = 1 << start_city, current_city = start_city;

        while (current_city != -1) {
            path.push_back(current_city);
            int next_city = parent[mask][current_city];
            mask |= (1 << next_city);
            current_city = next_city;
        }

        return path;
    };

    vector<int> min_cost_path = reconstructPath(parent_cost);
    vector<int> min_time_path = reconstructPath(parent_time);

    return {min_cost_path, min_time_path};
}

// Function to draw the flight paths using Dynamic Programming
void drawFlightDPPath(cv::Mat& base_img, const vector<cv::Point>& city_coords, const vector<int>& path, const string& output_file) {
    cv::Mat img = base_img.clone();

    // Draw cities and indices along the path
    for (size_t i = 0; i < path.size(); ++i) {
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
double calculateTotalMetric(const vector<int>& path, const unordered_map<string, vector<double>>& flight_data, const vector<string>& city_names, bool cost_metric = true) {
    double total_metric = 0.0;

    for (size_t i = 0; i < path.size() - 1; ++i) {
        string key = city_names[path[i]] + "_" + city_names[path[i + 1]];
        string reverse_key = city_names[path[i + 1]] + "_" + city_names[path[i]];

        int index = cost_metric ? 0 : 1; // 0 for cost, 1 for time

        if (flight_data.find(key) != flight_data.end()) {
            total_metric += flight_data.at(key)[index];
        } else if (flight_data.find(reverse_key) != flight_data.end()) {
            total_metric += flight_data.at(reverse_key)[index];
        } else {
            cerr << "Warning: No flight data found for route " << key << " or " << reverse_key << endl;
        }
    }

    return total_metric;
}

// Function to write the city order and metrics to a CSV file
void writeCityOrderToCSV(const string& filename, const vector<int>& path, const vector<string>& city_names, const vector<cv::Point>& city_coords, int total_cost, int total_time) {
    ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << "Index,City Name,Coordinate (x,y)" << endl;

        for (size_t new_index = 0; new_index < path.size(); ++new_index) {
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

    auto flight_data = loadFlightData("./Dataset/Dataset_Flight.csv", city_names);
    cv::Mat base_img = cv::imread("./Image/Europe.png");
    if (flight_data.empty()) return -1;

    // Display available cities
    cout << "Available cities:" << endl;
    for (size_t i = 0; i < city_names.size(); ++i) {
        cout << i << ": " << city_names[i] << " (" << city_coords[i].x << ", " << city_coords[i].y << ")" << endl;
    }

    // Get the starting city from the user
    int start_city;
    cout << "Enter the index of the starting city: ";
    cin >> start_city;
    if (start_city < 0 || start_city >= city_names.size()) {
        cerr << "Error: Invalid starting city index." << endl;
        return -1;
    }

    // Get paths for different scenarios
    auto [min_cost_path_greedy, min_time_path_greedy] = flightGreedyPath(city_coords, flight_data, city_names, start_city);
    auto [min_cost_path_dc, min_time_path_dc] = flightDCPath(city_coords, flight_data, city_names, start_city);
    auto [min_cost_path_dp, min_time_path_dp] = flightDPPath(city_coords, flight_data, city_names, start_city);

    // Calculate metrics for Greedy paths
    double total_cost_greedy_cost_path = calculateTotalMetric(min_cost_path_greedy, flight_data, city_names, true);
    double total_time_greedy_cost_path = calculateTotalMetric(min_cost_path_greedy, flight_data, city_names, false);

    double total_cost_greedy_time_path = calculateTotalMetric(min_time_path_greedy, flight_data, city_names, true);
    double total_time_greedy_time_path = calculateTotalMetric(min_time_path_greedy, flight_data, city_names, false);

    // Calculate metrics for Divide and Conquer paths
    double total_cost_dc_cost_path = calculateTotalMetric(min_cost_path_dc, flight_data, city_names, true);
    double total_time_dc_cost_path = calculateTotalMetric(min_cost_path_dc, flight_data, city_names, false);

    double total_cost_dc_time_path = calculateTotalMetric(min_time_path_dc, flight_data, city_names, true);
    double total_time_dc_time_path = calculateTotalMetric(min_time_path_dc, flight_data, city_names, false);

    // Calculate metrics for Dynamic Programming paths
    double total_cost_dp_cost_path = calculateTotalMetric(min_cost_path_dp, flight_data, city_names, true);
    double total_time_dp_cost_path = calculateTotalMetric(min_cost_path_dp, flight_data, city_names, false);

    double total_cost_dp_time_path = calculateTotalMetric(min_time_path_dp, flight_data, city_names, true);
    double total_time_dp_time_path = calculateTotalMetric(min_time_path_dp, flight_data, city_names, false);

    // Draw and save paths for Greedy
    drawFlightGreedyPath(base_img, city_coords, min_cost_path_greedy, "./Image/E3Flight_Greedy_MinCost_Path.png");
    drawFlightGreedyPath(base_img, city_coords, min_time_path_greedy, "./Image/E3Flight_Greedy_MinTime_Path.png");

    // Draw and save paths for Divide and Conquer
    drawFlightDCPath(base_img, city_coords, min_cost_path_dc, "./Image/E3Flight_DC_MinCost_Path.png");
    drawFlightDCPath(base_img, city_coords, min_time_path_dc, "./Image/E3Flight_DC_MinTime_Path.png");

    // Draw and save paths for DP
    drawFlightDPPath(base_img, city_coords, min_cost_path_dp, "./Image/E3Flight_DP_MinCost_Path.png");
    drawFlightDPPath(base_img, city_coords, min_time_path_dp, "./Image/E3Flight_DP_MinTime_Path.png");

    // Write Greedy paths to CSV
    writeCityOrderToCSV("./Table/E3Flight_Greedy_MinCost_Path.csv", min_cost_path_greedy, city_names, city_coords, total_cost_greedy_cost_path, total_time_greedy_cost_path);
    writeCityOrderToCSV("./Table/E3Flight_Greedy_MinTime_Path.csv", min_time_path_greedy, city_names, city_coords, total_cost_greedy_time_path, total_time_greedy_time_path);

    // Write Divide and Conquer paths to CSV
    writeCityOrderToCSV("./Table/E3Flight_DC_MinCost_Path.csv", min_cost_path_dc, city_names, city_coords, total_cost_dc_cost_path, total_time_dc_cost_path);
    writeCityOrderToCSV("./Table/E3Flight_DC_MinTime_Path.csv", min_time_path_dc, city_names, city_coords, total_cost_dc_time_path, total_time_dc_time_path);

    // Write DP paths to CSV
    writeCityOrderToCSV("./Table/E3Flight_DP_MinCost_Path.csv", min_cost_path_dp, city_names, city_coords, total_cost_dp_cost_path, total_time_dp_cost_path);
    writeCityOrderToCSV("./Table/E3Flight_DP_MinTime_Path.csv", min_time_path_dp, city_names, city_coords, total_cost_dp_time_path, total_time_dp_time_path);

    cout << "Greedy Min Cost Path: Total Cost: $" << total_cost_greedy_cost_path << ", Total Time: " << total_time_greedy_cost_path << " mins" << endl;
    cout << "Greedy Min Time Path: Total Cost: $" << total_cost_greedy_time_path << ", Total Time: " << total_time_greedy_time_path << " mins" << endl;

    cout << "DC Min Cost Path: Total Cost: $" << total_cost_dc_cost_path << ", Total Time: " << total_time_dc_cost_path << " mins" << endl;
    cout << "DC Min Time Path: Total Cost: $" << total_cost_dc_time_path << ", Total Time: " << total_time_dc_time_path << " mins" << endl;

    cout << "DP Min Cost Path: Total Cost: $" << total_cost_dp_cost_path << ", Total Time: " << total_time_dp_cost_path << " mins" << endl;
    cout << "DP Min Time Path: Total Cost: $" << total_cost_dp_time_path << ", Total Time: " << total_time_dp_time_path << " mins" << endl;

    // Open the saved images using the system's default image viewer
#ifdef _WIN32
    system("start ./Image/E3Flight_Greedy_MinCost_Path.png");
    system("start ./Image/E3Flight_Greedy_MinTime_Path.png");
    system("start ./Image/E3Flight_DC_MinCost_Path.png");
    system("start ./Image/E3Flight_DC_MinTime_Path.png");
    system("start ./Image/E3Flight_DP_MinCost_Path.png");
    system("start ./Image/E3Flight_DP_MinTime_Path.png");
#elif __APPLE__
    system("open ./Image/E3Flight_Greedy_MinCost_Path.png");
    system("open ./Image/E3Flight_Greedy_MinTime_Path.png");
    system("open ./Image/E3Flight_DC_MinCost_Path.png");
    system("open ./Image/E3Flight_DC_MinTime_Path.png");
    system("open ./Image/E3Flight_DP_MinCost_Path.png");
    system("open ./Image/E3Flight_DP_MinTime_Path.png");
#elif __linux__
    system("xdg-open ./Image/E3Flight_Greedy_MinCost_Path.png");
    system("xdg-open ./Image/E3Flight_Greedy_MinTime_Path.png");
    system("xdg-open ./Image/E3Flight_DC_MinCost_Path.png");
    system("xdg-open ./Image/E3Flight_DC_MinTime_Path.png");
    system("xdg-open ./Image/E3Flight_DP_MinCost_Path.png");
    system("xdg-open ./Image/E3Flight_DP_MinTime_Path.png");
#endif

    return 0;
}