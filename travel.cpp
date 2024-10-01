#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

const int INF = 1e9;

struct City {
    cv::Point point;
    int index;
};

// Function to load city coordinates from the CSV file
vector<cv::Point> loadCityCoordinates(const string& coor_file)
{
    vector<cv::Point> city_list;
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

        city_list.push_back(cv::Point(x, y));
    }

    infile.close();
    return city_list;
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


// Greedy algorithm to find the shortest path based on distance and return the path indices
vector<int> shortestGreedyPath(const vector<cv::Point>& city_list)
{
    int n = city_list.size();
    vector<bool> visited(n, false);  // Track visited cities
    vector<int> shortest_greedy_index;  // Store the path of city indices
    int current_city = 0;  // Start at the first city
    shortest_greedy_index.push_back(current_city);
    visited[current_city] = true;

    for (int i = 1; i < n; i++)
    {
        double min_distance = numeric_limits<double>::max();
        int next_city = -1;

        for (int j = 0; j < n; j++)
        {
            if (!visited[j])
            {
                double dist = calculateDistance(city_list[current_city], city_list[j]);
                if (dist < min_distance)
                {
                    min_distance = dist;
                    next_city = j;
                }
            }
        }

        if (next_city != -1)
        {
            shortest_greedy_index.push_back(next_city);
            visited[next_city] = true;
            current_city = next_city;
        }
    }

    return shortest_greedy_index;
}

// Function to draw the Shortest Path using Greedy Algorithm (plot cities, draw paths, and annotate them with their greedy index)
void drawShortestGreedyPath(cv::Mat& img, const vector<cv::Point>& city_list, const vector<int>& path)
{
    // Plot the cities and annotate them with their new greedy index
    for (size_t i = 0; i < path.size(); i++)
    {
        int greedy_index = path[i];
        // Plot city as a red circle
        cv::circle(img, city_list[greedy_index], 10, cv::Scalar(0, 0, 255), -1);

        // Annotate the city with its new greedy index
        cv::putText(img, to_string(i), city_list[greedy_index], cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw the greedy path with green lines based on the calculated path
    for (size_t i = 0; i < path.size() - 1; i++)
    {
        cv::line(img, city_list[path[i]], city_list[path[i + 1]], cv::Scalar(30, 150, 80), 4);
    }

    // Close the greedy path (draw line from last city to first city in path)
    if (!path.empty())
    {
        cv::line(img, city_list[path[path.size() - 1]], city_list[path[0]], cv::Scalar(30, 150, 80), 4);
    }
}


// Function to implement the Divide-and-Conquer algorithm to find the closest pair
pair<int, int> divideAndConquer(vector<City>& cities, int left, int right, double& min_dist)
{
    // Base case: If there are 3 or fewer cities, use brute-force
    if (right - left <= 1)
    {
        double dist = calculateDistance(cities[left].point, cities[right].point);
        min_dist = dist;
        return { cities[left].index, cities[right].index };
    }

    // Find the middle point
    int mid = left + (right - left) / 2;
    cv::Point mid_point = cities[mid].point;

    // Recursively find the closest pair on both halves
    double left_dist = numeric_limits<double>::max(), right_dist = numeric_limits<double>::max();
    pair<int, int> left_pair = divideAndConquer(cities, left, mid, left_dist);
    pair<int, int> right_pair = divideAndConquer(cities, mid + 1, right, right_dist);

    // Determine the smallest distance from the left and right halves
    min_dist = min(left_dist, right_dist);
    pair<int, int> best_pair = (min_dist == left_dist) ? left_pair : right_pair;

    return best_pair;
}

// Function to implement the Divide-and-Conquer strategy to find the shortest path and return the path indices
vector<int> shortestDCPath(const vector<cv::Point>& city_list)
{
    // Convert cv::Point to City struct for DC algorithm
    vector<City> city_list_dc;
    for (size_t i = 0; i < city_list.size(); i++)
    {
        city_list_dc.push_back({city_list[i], static_cast<int>(i)});
    }

    // Sort cities by x-coordinates for the divide and conquer approach
    sort(city_list_dc.begin(), city_list_dc.end(), [](const City& a, const City& b) { return a.point.x < b.point.x; });

    double min_dist = numeric_limits<double>::max();
    pair<int, int> closest_pair = divideAndConquer(city_list_dc, 0, city_list_dc.size() - 1, min_dist);

    vector<int> dc_path;
    dc_path.push_back(closest_pair.first);
    dc_path.push_back(closest_pair.second);

    // Sort the cities based on y-coordinates for visualizing the path
    sort(city_list_dc.begin(), city_list_dc.end(), [](const City& a, const City& b) { return a.point.y < b.point.y; });

    for (const auto& city : city_list_dc)
    {
        if (city.index != closest_pair.first && city.index != closest_pair.second)
            dc_path.push_back(city.index);
    }

    return dc_path;
}

// Function to draw the Shortest Path using Divide-and-Conquer and save as an image
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

// Function to compute shortest path using Dynamic Programming (TSP)
vector<int> shortestDPPath(const vector<cv::Point>& city_list)
{
    int n = city_list.size();
    vector<vector<double>> dist(n, vector<double>(n, 0));

    // Fill the distance matrix
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            dist[i][j] = calculateDistance(city_list[i], city_list[j]);
        }
    }

    // DP table where dp[mask][i] is the minimum distance to visit all cities in 'mask' and end at city i
    vector<vector<double>> dp(1 << n, vector<double>(n, INF));
    vector<vector<int>> parent(1 << n, vector<int>(n, -1)); // To reconstruct the path

    // Start from city 0
    dp[1][0] = 0;

    // Iterate over all subsets of cities
    for (int mask = 1; mask < (1 << n); mask++)
    {
        for (int u = 0; u < n; u++)
        {
            if (!(mask & (1 << u))) continue; // Skip if u is not in mask

            // Try to go to the next city v
            for (int v = 0; v < n; v++)
            {
                if (mask & (1 << v)) continue; // Skip if v is already visited

                int next_mask = mask | (1 << v);
                double new_dist = dp[mask][u] + dist[u][v];

                if (new_dist < dp[next_mask][v])
                {
                    dp[next_mask][v] = new_dist;
                    parent[next_mask][v] = u; // Record the path
                }
            }
        }
    }

    // Reconstruct the path (find minimum distance to return to city 0)
    double min_cost = INF;
    int last_city = -1;
    int final_mask = (1 << n) - 1;

    for (int i = 1; i < n; i++)
    {
        double final_cost = dp[final_mask][i] + dist[i][0];
        if (final_cost < min_cost)
        {
            min_cost = final_cost;
            last_city = i;
        }
    }

    // Reconstruct the shortest path
    vector<int> path;
    int current_city = last_city;
    int current_mask = final_mask;

    while (current_city != -1)
    {
        path.push_back(current_city);
        int temp = parent[current_mask][current_city];
        current_mask ^= (1 << current_city);
        current_city = temp;
    }

    path.push_back(0); // Add the starting city to complete the cycle
    reverse(path.begin(), path.end()); // Reverse to get the correct order

    return path;
}

// Function to draw the Shortest Path using Dynamic Programming (TSP)
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

int main()
{
    // Load the city coordinates from CSV file
    string coor_file = "./Dataset/Dataset_Coordinate.csv";
    vector<cv::Point> city_list = loadCityCoordinates(coor_file);
    if (city_list.empty()) return -1;  // Exit if failed to load city coordinates

    // Load the image (replace with your image path)
    string image_file = "./Image/Europe.png";
    cv::Mat img_initial = cv::imread(image_file);  // For initial path
    cv::Mat img_greedy = cv::imread(image_file);   // For greedy path
    cv::Mat img_dc = cv::imread(image_file);       // For divide-and-conquer path
    cv::Mat img_dp = cv::imread(image_file);       // For dynamic programming path

    if (img_initial.empty() || img_greedy.empty() || img_dc.empty())
    {
        cerr << "Error: Unable to load image." << endl;
        return -1;
    }

    // 1. Draw the initial path on img_initial
    drawInitialPath(img_initial, city_list);

    // 2. Apply the greedy algorithm to find the shortest path
    vector<int> shortest_greedy_index = shortestGreedyPath(city_list);

    // 2. Draw the greedy path on img_greedy and plot cities with their greedy index
    drawShortestGreedyPath(img_greedy, city_list, shortest_greedy_index);

    // 3. Apply the divide-and-conquer algorithm to find the shortest path
    vector<int> shortest_dc_index = shortestDCPath(city_list);

    // 3. Draw the divide-and-conquer path on img_dc and plot cities with their index
    drawShortestDCPath(img_dc, city_list, shortest_dc_index);

    // 4. Apply the dynamic programming algorithm to find the shortest path
    vector<int> shortest_dp_index = shortestDPPath(city_list);

    // 4. Draw the dynamic programming path on img_dp and plot cities with their DP index
    drawShortestDPPath(img_dp, city_list, shortest_dp_index);
    
    string initial_path_file = "./Image/Initial_Path.png";
    string shortest_greedy_path_file = "./Image/Shortest_Greedy_Path.png";
    string shortest_dc_path_file = "./Image/Shortest_DC_Path.png";
    string shortest_dp_path_file = "./Image/Shortest_DP_Path.png";

    cv::imwrite(initial_path_file, img_initial);
    cv::imwrite(shortest_greedy_path_file, img_greedy);
    cv::imwrite(shortest_dc_path_file, img_dc);
    cv::imwrite(shortest_dp_path_file, img_dp);

    cout << "Initial path saved to " << initial_path_file << endl;
    cout << "Shortest path (greedy) saved to " << shortest_greedy_path_file << endl;
    cout << "Shortest path (divide-and-conquer) saved to " << shortest_dc_path_file << endl;
    cout << "Shortest path (dynamic programming) saved to " << shortest_dp_path_file << endl;

    // Open the saved images using the system's default image viewer
#ifdef _WIN32
    system(("start " + initial_path_file).c_str());
    system(("start " + shortest_greedy_path_file).c_str());
    system(("start " + shortest_dc_path_file).c_str());
    system(("start " + shortest_dp_path_file).c_str());
#elif __APPLE__
    system(("open " + initial_path_file).c_str());
    system(("open " + shortest_greedy_path_file).c_str());
    system(("open " + shortest_dc_path_file).c_str());
    system(("open " + shortest_dp_path_file).c_str());
#elif __linux__
    system(("xdg-open " + initial_path_file).c_str());
    system(("xdg-open " + shortest_greedy_path_file).c_str());
    system(("xdg-open " + shortest_dc_path_file).c_str());
    system(("xdg-open " + shortest_dp_path_file).c_str());
#endif

    return 0;
}
