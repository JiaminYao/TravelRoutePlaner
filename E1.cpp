#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <unordered_set>

using namespace std;

const int INF = 1e9;

struct City {
    cv::Point point;
    int index;
    string name;
};

// Function to load city coordinates from the CSV file
vector<City> loadCityCoordinates(const string& coor_file)
{
    vector<City> city_list;
    ifstream infile(coor_file);

    if (!infile.is_open())
    {
        cerr << "Error: Could not open the coordinates file." << endl;
        return city_list;
    }

    string line;
    getline(infile, line);  // Skip the header line if it exists

    // Read the CSV file line by line
    while (getline(infile, line))
    {
        stringstream ss(line);
        string city_name, x_str, y_str;
        int x, y;

        // Extract city name, x-coordinate, and y-coordinate
        getline(ss, city_name, ',');
        getline(ss, x_str, ',');
        getline(ss, y_str, ',');

        try {
            x = stoi(x_str);
            y = stoi(y_str);
        } catch (invalid_argument&)
        {
            cerr << "Error: Invalid numeric value in coordinates for city " << city_name << endl;
            continue;
        }

        city_list.push_back({cv::Point(x, y), static_cast<int>(city_list.size()), city_name});
    }

    infile.close();
    return city_list;
}


// Function to write city names and their order based on the calculated path to a CSV file
void writeCityOrderToCSV(const string& filename, const vector<City>& city_list, const vector<int>& path)
{
    ofstream outfile(filename);
    if (outfile.is_open()) {

    // Write header
    outfile << "Index,City Name,Coordinate (x,y)" << endl;

    // Write city information based on the path
    for (size_t new_index = 0; new_index < path.size(); ++new_index)
    {
        int original_index = path[new_index];
        outfile << new_index << "," << city_list[original_index].name 
                << ",(" << city_list[original_index].point.x 
                << "," << city_list[original_index].point.y << ")" << endl;
    }
    outfile.close();
    } else {
        cerr << "Error: Could not open file " << filename << " for writing." << endl;
    }
}

// Function to calculate the distance between two points
double calculateDistance(const cv::Point& a, const cv::Point& b)
{
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}


// Function to draw Initial Path without order (plot cities, draw paths, and annotate them with their index)
void drawInitialPath(cv::Mat& img, const vector<cv::Point>& city_list)
{
    // Plot the cities, draw paths, and annotate with index
    for (size_t i = 0; i < city_list.size(); i++)
    {
        // Plot city as a red circle
        cv::circle(img, city_list[i], 10, cv::Scalar(0, 0, 255), -1);

        // Annotate the city with its index
        cv::putText(img, to_string(i), city_list[i], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);

        // Draw line to the next city (if there is one)
        if (i < city_list.size() - 1)
        {
            cv::line(img, city_list[i], city_list[i + 1], cv::Scalar(255, 0, 0), 4);
        }
    }

    // Close the path (draw line from last city to first city)
    if (!city_list.empty())
    {
        cv::line(img, city_list[city_list.size() - 1], city_list[0], cv::Scalar(255, 0, 0), 4);
    }
}


// Function to implement Greedy algorithm to find the shortest path
vector<int> shortestGreedyPath(const vector<cv::Point>& city_list) {
    int n = city_list.size();
    vector<bool> visited(n, false);  // Track visited cities
    vector<int> shortest_greedy_index;  // Store the path of city indices
    int current_city = 0;  // Start at the first city
    shortest_greedy_index.push_back(current_city);
    visited[current_city] = true;

    for (int i = 1; i < n; i++) {
        double min_distance = numeric_limits<double>::max();
        int next_city = -1;

        for (int j = 0; j < n; j++) {
            if (!visited[j]) {
                double dist = calculateDistance(city_list[current_city], city_list[j]);
                if (dist < min_distance) {
                    min_distance = dist;
                    next_city = j;
                }
            }
        }

        if (next_city != -1) {
            shortest_greedy_index.push_back(next_city);
            visited[next_city] = true;
            current_city = next_city;
        }
    }

    return shortest_greedy_index;
}

// Function to draw the Shortest Path using Greedy Algorithm (plot cities, draw paths, and annotate them with their greedy index)
void drawShortestGreedyPath(cv::Mat& img, const vector<cv::Point>& city_coords, const vector<int>& path)
{
    // Plot the cities and annotate them with their new greedy index
    for (size_t i = 0; i < path.size(); i++)
    {
        int greedy_index = path[i];
        // Plot city as a red circle
        cv::circle(img, city_coords[greedy_index], 10, cv::Scalar(0, 0, 255), -1);

        // Annotate the city with its new greedy index
        cv::putText(img, to_string(i), city_coords[greedy_index], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw the greedy path with green lines based on the calculated path
    for (size_t i = 0; i < path.size() - 1; i++)
    {
        cv::line(img, city_coords[path[i]], city_coords[path[i + 1]], cv::Scalar(30, 150, 80), 4);
    }

    // Close the greedy path (draw line from last city to first city in path)
    if (!path.empty())
    {
        cv::line(img, city_coords[path[path.size() - 1]], city_coords[path[0]], cv::Scalar(30, 150, 80), 4);
    }
}


// Function to implement the Divide-and-Conquer algorithm to find the shortest pairs
void divideAndConquerHelper(const vector<City>& cities, vector<int>& path) {
    if (cities.empty()) return;

    if (cities.size() == 1) {
        // Add the single city's index to the path
        path.push_back(cities[0].index);
        return;
    }

    size_t mid = cities.size() / 2;

    // Divide the cities into two halves
    vector<City> left_cities(cities.begin(), cities.begin() + mid);
    vector<City> right_cities(cities.begin() + mid, cities.end());

    // Process left and right halves
    divideAndConquerHelper(left_cities, path);
    divideAndConquerHelper(right_cities, path);
}


// Main function for Divide-and-Conquer
vector<int> shortestDCPath(const vector<City>& city_list, int start_city_index) {

    // Validate the start city index
    if (start_city_index < 0 || start_city_index >= city_list.size()) {
        cerr << "Error: Invalid start city index!" << endl;
        return {};
    }

    // Identify the starting city
    City start_city = city_list[start_city_index];

    // Create a new list excluding the starting city
    vector<City> remaining_cities;
    for (int i = 0; i < city_list.size(); ++i) {
        if (i != start_city_index) {
            remaining_cities.push_back(city_list[i]);
        }
    }

    // Sort the remaining cities by y-coordinate
    sort(remaining_cities.begin(), remaining_cities.end(), [](const City& a, const City& b) {
        return a.point.y < b.point.y;
    });

    // Recursive Divide-and-Conquer
    function<void(const vector<City>&, vector<int>&)> divideAndConquerHelper = [&](const vector<City>& cities, vector<int>& path) {
        if (cities.size() <= 1) {
            if (!cities.empty()) {
                path.push_back(cities[0].index); // Add the single city to the path
            }
            return;
        }

        size_t mid = cities.size() / 2;

        // Divide the cities into two halves
        vector<City> left_cities(cities.begin(), cities.begin() + mid);
        vector<City> right_cities(cities.begin() + mid, cities.end());

        // Recur on both halves
        divideAndConquerHelper(left_cities, path);
        divideAndConquerHelper(right_cities, path);
    };

    // Start processing
    vector<int> path;
    path.push_back(start_city.index); // Add the starting city first
    divideAndConquerHelper(remaining_cities, path);

    return path;
}

// Function to draw the Shortest Path using Divide-and-Conquer Algorithm
void drawShortestDCPath(cv::Mat& img, const vector<cv::Point>& city_list, const vector<int>& dc_path)
{
    // Plot the cities and annotate them with their index
    for (size_t i = 0; i < dc_path.size(); i++)
    {
        int index = dc_path[i];
        // Plot city as a red circle
        cv::circle(img, city_list[index], 10, cv::Scalar(0, 0, 255), -1);

        // Annotate the city with its new DC index
        cv::putText(img, to_string(i), city_list[index], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw the DC path with blue lines based on the calculated path
    for (size_t i = 0; i < dc_path.size() - 1; i++)
    {
        cv::line(img, city_list[dc_path[i]], city_list[dc_path[i + 1]], cv::Scalar(170, 30, 220), 4);
    }

    // Close the DC path (draw line from last city to first city in path)
    if (!dc_path.empty())
    {
        cv::line(img, city_list[dc_path[dc_path.size() - 1]], city_list[dc_path[0]], cv::Scalar(170, 30, 220), 4);
    }
}

// Function to implement the Dynamic Programming algorithm to find the shortest path
vector<int> shortestDPPath(const vector<cv::Point>& city_list) {
    int n = city_list.size();
    vector<vector<double>> dist(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dist[i][j] = calculateDistance(city_list[i], city_list[j]);
        }
    }

    vector<vector<double>> dp(1 << n, vector<double>(n, INF));
    vector<vector<int>> parent(1 << n, vector<int>(n, -1));
    dp[1][0] = 0;

    for (int mask = 1; mask < (1 << n); mask++) {
        for (int u = 0; u < n; u++) {
            if (!(mask & (1 << u))) continue;
            for (int v = 0; v < n; v++) {
                if (mask & (1 << v)) continue;
                int next_mask = mask | (1 << v);
                double new_dist = dp[mask][u] + dist[u][v];
                if (new_dist < dp[next_mask][v]) {
                    dp[next_mask][v] = new_dist;
                    parent[next_mask][v] = u;
                }
            }
        }
    }

    vector<int> path;
    double min_cost = INF;
    int last_city = -1;
    int final_mask = (1 << n) - 1;

    for (int i = 1; i < n; i++) {
        double cost = dp[final_mask][i] + dist[i][0];
        if (cost < min_cost) {
            min_cost = cost;
            last_city = i;
        }
    }

    int current_city = last_city;
    int current_mask = final_mask;

    while (current_city != -1) {
        path.push_back(current_city);
        int temp = parent[current_mask][current_city];
        current_mask ^= (1 << current_city);
        current_city = temp;
    }

    reverse(path.begin(), path.end());
    return path;
}

// Function to draw the Shortest Path using Dynamic Programming Algorithm
void drawShortestDPPath(cv::Mat& img, const vector<cv::Point>& city_list, const vector<int>& dp_path)
{
    // Plot the cities and annotate them with their index
    for (size_t i = 0; i < dp_path.size(); i++)
    {
        int index = dp_path[i];
        // Plot city as a red circle
        cv::circle(img, city_list[index], 10, cv::Scalar(0, 0, 255), -1);

        // Annotate the city with its DP index
        cv::putText(img, to_string(i), city_list[index], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw the DP path with yellow lines based on the calculated path
    for (size_t i = 0; i < dp_path.size() - 1; i++)
    {
        cv::line(img, city_list[dp_path[i]], city_list[dp_path[i + 1]], cv::Scalar(0, 255, 255), 4);
    }

    // Close the DP path (draw line from last city to first city in path)
    if (!dp_path.empty())
    {
        cv::line(img, city_list[dp_path[dp_path.size() - 1]], city_list[dp_path[0]], cv::Scalar(0, 255, 255), 4);
    }
}


// Function to calculate the total distance of a path
double calculateTotalDistance(const vector<cv::Point>& city_list, const vector<int>& path)
{
    double total_distance = 0.0;
    
    // Sum the distances between consecutive cities
    for (size_t i = 0; i < path.size() - 1; i++)
    {
        total_distance += calculateDistance(city_list[path[i]], city_list[path[i + 1]]);
    }

    // Add the distance from the last city back to the starting city
    total_distance += calculateDistance(city_list[path.back()], city_list[path[0]]);

    return total_distance;
}


int main()
{
    // Load the city coordinates from CSV file
    string coor_file = "./Dataset/Dataset_Coordinate.csv";
    vector<City> city_coords = loadCityCoordinates(coor_file);
    if (city_coords.empty()) {
        cerr << "Error: City coordinates file is empty or not loaded." << endl;
        return -1;
    }

    // Load the image (replace with your image path)
    string image_file = "./Image/Europe.png";
    cv::Mat img_initial = cv::imread(image_file);  // For initial path
    cv::Mat img_greedy = cv::imread(image_file);   // For greedy path
    cv::Mat img_dc = cv::imread(image_file);       // For divide-and-conquer path
    cv::Mat img_dp = cv::imread(image_file);       // For dynamic programming path

    if (img_initial.empty() || img_greedy.empty() || img_dc.empty() || img_dp.empty())
    {
        cerr << "Error: Unable to load image." << endl;
        return -1;
    }

    // Display available cities
    cout << "Available cities:" << endl;
    for (size_t i = 0; i < city_coords.size(); ++i) {
        cout << i << ": " << city_coords[i].name << " (" << city_coords[i].point.x << ", " << city_coords[i].point.y << ")" << endl;
    }

    int start_city_index;
    cout << "Enter the index of the starting city: ";
    cin >> start_city_index;

    if (start_city_index < 0 || start_city_index >= city_coords.size()) {
        cerr << "Error: Invalid city index." << endl;
        return -1;
    }

    // 1. **Initial Path** (Just reorder to show the starting city first)
    vector<cv::Point> city_coordinates;
    for (const auto& city : city_coords) {
        city_coordinates.push_back(city.point);
    }
    vector<int> initial_path;
    initial_path.push_back(start_city_index);
    for (size_t i = 0; i < city_coords.size(); ++i) {
        if (static_cast<int>(i) != start_city_index) {
            initial_path.push_back(i);
        }
    }
    drawInitialPath(img_initial, city_coordinates);

    // 2. **Greedy Algorithm**
    vector<cv::Point> greedy_input;
    for (int index : initial_path) {
        greedy_input.push_back(city_coordinates[index]);
    }
    vector<int> greedy_path = shortestGreedyPath(greedy_input);
    vector<int> mapped_greedy_path;
    for (int idx : greedy_path) {
        mapped_greedy_path.push_back(initial_path[idx]);
    }
    drawShortestGreedyPath(img_greedy, city_coordinates, mapped_greedy_path);
    double greedy_total_distance = calculateTotalDistance(city_coordinates, mapped_greedy_path);
    cout << "Total distance (Greedy): " << greedy_total_distance << " units = "
         << greedy_total_distance * 0.9 << " km" << endl;

    // 3. **Divide-and-Conquer Algorithm**
    vector<int> dc_path = shortestDCPath(city_coords, start_city_index);
    drawShortestDCPath(img_dc, city_coordinates, dc_path);
    double dc_total_distance = calculateTotalDistance(city_coordinates, dc_path);
    cout << "Total distance (Divide-and-Conquer): " << dc_total_distance << " units = "
         << dc_total_distance * 0.9 << " km" << endl;

    // 4. **Dynamic Programming Algorithm** (Fix for start city)
    vector<cv::Point> dp_input;
    for (int index : initial_path) {
        dp_input.push_back(city_coordinates[index]);
    }
    vector<int> dp_path = shortestDPPath(dp_input);
    vector<int> mapped_dp_path;
    for (int idx : dp_path) {
        mapped_dp_path.push_back(initial_path[idx]);
    }
    drawShortestDPPath(img_dp, city_coordinates, mapped_dp_path);
    double dp_total_distance = calculateTotalDistance(city_coordinates, mapped_dp_path);
    cout << "Total distance (Dynamic Programming): " << dp_total_distance << " units = "
         << dp_total_distance * 0.9 << " km" << endl;

    // Save images for paths
    cv::imwrite("./Image/Initial_Path.png", img_initial);
    cv::imwrite("./Image/E1Helicopter_Greedy_Path.png", img_greedy);
    cv::imwrite("./Image/E1Helicopter_DC_Path.png", img_dc);
    cv::imwrite("./Image/E1Helicopter_DP_Path.png", img_dp);

    // Save paths to CSV
    writeCityOrderToCSV("./Table/E1Helicopter_Greedy_Order.csv", city_coords, mapped_greedy_path);
    writeCityOrderToCSV("./Table/E1Helicopter_DC_Order.csv", city_coords, dc_path);
    writeCityOrderToCSV("./Table/E1Helicopter_DP_Order.csv", city_coords, mapped_dp_path);
    
    // Open the saved images using the system's default image viewer
#ifdef _WIN32
    system("start ./Image/Initial_Path.png");
    system("start ./Image/E1Helicopter_Greedy_Path.png");
    system("start ./Image/E1Helicopter_DC_Path.png");
    system("start ./Image/E1Helicopter_DP_Path.png");
#elif __APPLE__
    system("open ./Image/Initial_Path.png");
    system("open ./Image/E1Helicopter_Greedy_Path.png");
    system("open ./Image/E1Helicopter_DC_Path.png");
    system("open ./Image/E1Helicopter_DP_Path.png");
#elif __linux__
    system("xdg-open ./Image/Initial_Path.png");
    system("xdg-open ./Image/E1Helicopter_Greedy_Path.png");
    system("xdg-open ./Image/E1Helicopter_DC_Path.png");
    system("xdg-open ./Image/E1Helicopter_DP_Path.png");
#endif

    return 0;
}
