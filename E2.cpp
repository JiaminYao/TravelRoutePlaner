#include <opencv2/opencv.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <cmath>
#include <unordered_map>

using namespace std;

struct RouteInfo {
    string city1, city2;
    double regular_cost, regular_time;
    double highway_cost, highway_time;
};

struct CityVisit {
    string city;
    double cost;
    double time;
    cv::Point point;
};

unordered_map<string, cv::Point> loadCityCoordinates(const string& filename) {
    unordered_map<string, cv::Point> cityCoordinates;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error: Could not open the file " << filename << endl;
        return cityCoordinates; // Return empty map if file cannot be opened
    }

    string line;
    getline(file, line); // Skip the header line

    while (getline(file, line)) {
        istringstream ss(line);
        string city;
        int x, y;

        // Parse the city name and coordinates
        getline(ss, city, ',');
        ss >> x;
        ss.ignore(1, ','); // Skip the comma
        ss >> y;

        // Store the data in the unordered_map
        cityCoordinates[city] = cv::Point(x, y);
    }

    return cityCoordinates;
}

vector<RouteInfo> loadRoutesInfo(const string &filename) {
    vector<RouteInfo> routes;
    ifstream file(filename);
    string title;
    string line;

    if (!file.is_open()) {
        cerr << "Error: Could not open the file " << filename << endl;
        return routes; // Return empty vector if file cannot be opened
    }

    // Read the title line (header)
    getline(file, title);

    // Read each subsequent line
    while (getline(file, line)) {
        istringstream ss(line);
        RouteInfo route;

        // Read each field, separated by commas
        getline(ss, route.city1, ',');
        getline(ss, route.city2, ',');
        ss >> route.regular_cost;
        ss.ignore(); // Ignore the comma
        ss >> route.regular_time;
        ss.ignore(); // Ignore the comma
        ss >> route.highway_cost;
        ss.ignore(); // Ignore the comma
        ss >> route.highway_time;

        routes.push_back(route);
    }

    file.close(); // Close the file
    return routes;
}

// Function to check if a character is a space
bool isWhitespace(char c) {
    return isspace(static_cast<unsigned char>(c)); // Cast to unsigned char for safety
}

vector<string> loadCities(const string &filename) {
    vector<string> cities;
    ifstream file(filename);
    string title;
    getline(file, title);

    if (!file.is_open()) {
        cerr << "Error: Could not open the file " << filename << endl;
        return cities; // Return empty vector if file cannot be opened
    }

    string city;
    while (getline(file, city)) {
        // Remove any unwanted newline characters and trim whitespace
        city.erase(remove(city.begin(), city.end(), '\n'), city.end()); // Remove newlines
        city.erase(remove_if(city.begin(), city.end(), isWhitespace), city.end()); // Remove spaces
        cities.push_back(city);
    }

    file.close(); // Close the file
    return cities;
}

// Function to find the next city using a greedy approach
CityVisit findNextCity(const string& current_city, vector<RouteInfo>& roads, unordered_map<string, bool>& visited, unordered_map<string, cv::Point>& citycoordinates, bool isHighway) {
    CityVisit next_city_data = {"", numeric_limits<double>::max(), numeric_limits<double>::max()};
    //cout << roads[0].city1 << endl;
    //cout << current_city<< endl;
    for (const auto& road : roads) {
        
        if (road.city1 == current_city && !visited[road.city2] ) {
            double cost = isHighway ? road.highway_cost : road.regular_cost;
            double time = isHighway ? road.highway_time : road.regular_time;
            if (cost < next_city_data.cost || (cost == next_city_data.cost && time < next_city_data.time)) {
                next_city_data.city = road.city2;
                next_city_data.cost = cost;
                next_city_data.time = time;
                next_city_data.point = citycoordinates[road.city2];
            }
        } else if (road.city2 == current_city && !visited[road.city1]) {
            double cost = isHighway ? road.highway_cost : road.regular_cost;
            double time = isHighway ? road.highway_time : road.regular_time;

            if (cost < next_city_data.cost || (cost == next_city_data.cost && time < next_city_data.time)) {
                next_city_data.city = road.city1;
                next_city_data.cost = cost;
                next_city_data.time = time;
                next_city_data.point = citycoordinates[road.city1];
            }
        }
    }

    return next_city_data;
}

// Greedy algorithm to minimize total cost and time
vector<CityVisit> travelGreedy(const vector<string>& cities_to_visit, vector<RouteInfo>& roads, unordered_map<string, cv::Point>& citycoordinates, bool isHighway) {
    unordered_map<string, bool> visited;
    vector<CityVisit> travel_sequence;

    // Initialize visited cities map
    for (const auto& city : cities_to_visit) {
        visited[city] = false;
    }

    // Start from the first city in the cities_to_visit list
    string current_city = cities_to_visit[0];
    visited[current_city] = true;
    CityVisit c;
    c.city = current_city, c.cost = 0, c.time = 0;
    travel_sequence.push_back(c);

    for (size_t i = 1; i < cities_to_visit.size(); ++i) {
        CityVisit next_city_data = findNextCity(current_city, roads, visited, citycoordinates, isHighway);
        //cout << i << endl;
        if (next_city_data.city.empty()) {
            cout << "No path found to complete the journey." << endl;
            return vector<CityVisit>();
        }

        travel_sequence.push_back(next_city_data);
        visited[next_city_data.city] = true;
        current_city = next_city_data.city;
    }

    return travel_sequence;
}

// Function to print the travel sequence
void printTravelSequence(const vector<CityVisit>& sequence) {
    double total_cost = 0, total_time = 0;
    cout << "Travel Sequence:" << endl;
    for (const auto& step : sequence) {
        cout << "City: " << step.city << ", Cost: " << step.cost << ", Time: " << step.time << endl;
        total_cost += step.cost;
        total_time += step.time;
    }
    cout << "Total Cost: " << total_cost << ", Total Time: " << total_time << endl;
}

/* void writeCityOrderToCSV(const string& file_name, const vector<City>& city_list, const vector<int>& path)
{
    ofstream outfile(file_name);
    if (!outfile.is_open())
    {
        cerr << "Error: Could not open the file " << file_name << " for writing." << endl;
        return;
    }

    // Write header
    outfile << "Index,City Name,Coordinate (x,y)\n";

    // Write city information based on the path
    for (size_t i = 0; i < path.size(); i++)
    {
        int idx = path[i];
        outfile << i << "," << city_list[idx].name << ",(" << city_list[idx].point.x << "," << city_list[idx].point.y << ")\n";
    }

    outfile.close();
    cout << "City order saved to " << file_name << endl;
}
 */
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

int main() {
    //file names for dataset
    string filenameD = "Dataset/Dataset_Driving.csv";
    string filenameC = "Dataset/Cites.csv";
    string filenameCC = "Dataset/Dataset_Coordinate.csv";
    
    //load the data seperatly, initialize the vector that store the final path
    vector<string> cities = loadCities(filenameC);
    vector<RouteInfo> drivingInfo = loadRoutesInfo(filenameD);
    unordered_map<string, cv::Point> cityCoordinates = loadCityCoordinates(filenameCC);
    vector<CityVisit> visitListR;
    vector<CityVisit> visitListH;
    
    // Load the image (replace with your image path)
    string image_file = "Image/Europe.png";
    cv::Mat img_initial = cv::imread(image_file);  // For initial path
    cv::Mat img_greedy = cv::imread(image_file);   // For greedy path
    cv::Mat img_dc = cv::imread(image_file);       // For divide-and-conquer path
    cv::Mat img_dp = cv::imread(image_file);       // For dynamic programming path

    if (img_initial.empty() || img_greedy.empty() || img_dc.empty() || img_dp.empty())
    {
        cerr << "Error: Unable to load image." << endl;
        return -1;
    }
    cout << 7999999999 << endl;
    vector<cv::Point> coordinates;
    for (const auto& city : cityCoordinates) {
        coordinates.push_back(city.second); // Only push back the coordinates (cv::Point)
        cout << 7 << endl;
    }
    

    // Assuming drawInitialPath is defined and takes an image and a vector of Points
    drawInitialPath(img_initial, coordinates);

    // Run greedy algorithm for regular roads
    visitListR = travelGreedy(cities, drivingInfo, cityCoordinates, false);
    printTravelSequence(visitListR);

    // Greedy solution for both regular and highway
    //greedyPath(cities, roads, false, greedyPathReg);
    //greedyPath(cities, roads, true, greedyPathHwy);

    // Dynamic Programming solution for both regular and highway
    //dynamicProgrammingPath(cities, roads, false, dpPathReg);
    //dynamicProgrammingPath(cities, roads, true, dpPathHwy);

    // Divide and Conquer solution for both regular and highway
    //divideAndConquerPath(cities, roads, false, dcPathReg);
    //divideAndConquerPath(cities, roads, true, dcPathHwy);

    // Output Results (greedy, DP, DC paths for regular and highway as required)
    string initial_path_file = "Image/E2Initial_Path.png";
    string shortest_greedy_path_file = "./Image/E2Shortest_Greedy_Path.png";
    string shortest_dc_path_file = "./Image/E2hortest_DC_Path.png";
    string shortest_dp_path_file = "./Image/E2Shortest_DP_Path.png";

    cv::imwrite(initial_path_file, img_initial);

    cout << "Initial path saved to " << initial_path_file << endl;


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