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

// Define path information for each type
struct PathInfo {
    cv::Mat map;
    string outputname;
    string path_name;
    //string image_filename;
    string csv_filename;
};


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

struct PathResult {
    double cost;
    double time;
    vector<string> sequence; // Stores the travel sequence of cities
    vector<pair<double, double>> costTimeAtEachStep; // Stores cost and time for each step
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
CityVisit findNextCity(const string& current_city, vector<RouteInfo>& roads, unordered_map<string, bool>& visited, unordered_map<string, cv::Point>& citycoordinates, bool isHighway, bool isMinCost) {
    CityVisit next_city_data = {"", numeric_limits<double>::max(), numeric_limits<double>::max()};
    
    for (const auto& road : roads) {
        // Only consider unvisited cities in the direction from current_city
        if (road.city1 == current_city && !visited[road.city2]) {
            double cost = isHighway ? road.highway_cost : road.regular_cost;
            double time = isHighway ? road.highway_time : road.regular_time;
            
            if (isMinCost) {
                // Prioritize by minimum cost, with time as a tiebreaker
                if (cost < next_city_data.cost || (cost == next_city_data.cost && time < next_city_data.time)) {
                    next_city_data.city = road.city2;
                    next_city_data.cost = cost;
                    next_city_data.time = time;
                    next_city_data.point = citycoordinates[road.city2];
                }
            } else {
                // Prioritize by minimum time, with cost as a tiebreaker
                if (time < next_city_data.time || (time == next_city_data.time && cost < next_city_data.cost)) {
                    next_city_data.city = road.city2;
                    next_city_data.cost = cost;
                    next_city_data.time = time;
                    next_city_data.point = citycoordinates[road.city2];
                }
            }
        }
    }
    
    return next_city_data;
}


// Greedy algorithm to minimize total cost and time
vector<CityVisit> travelGreedy(const vector<string>& cities_to_visit, vector<RouteInfo>& roads, unordered_map<string, cv::Point>& citycoordinates, bool isHighway, bool isMinCost) {
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
    c.city = current_city, c.cost = 0, c.time = 0, c.point = citycoordinates[current_city];
    travel_sequence.push_back(c);

    for (size_t i = 1; i < cities_to_visit.size(); ++i) {
        CityVisit next_city_data = findNextCity(current_city, roads, visited, citycoordinates, isHighway, isMinCost);
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


vector<CityVisit> tspDP(const vector<string>& cities, const vector<RouteInfo>& drivingInfo, unordered_map<string, cv::Point>& citycoordinates, bool useHighway, bool minimizeCost) {
    int n = cities.size();
    int allVisited = (1 << n) - 1;

    unordered_map<string, int> cityIndex;
    for (int i = 0; i < n; ++i) {
        cityIndex[cities[i]] = i;
    }

    vector<vector<double>> costMatrix(n, vector<double>(n, numeric_limits<double>::max()));
    vector<vector<double>> timeMatrix(n, vector<double>(n, numeric_limits<double>::max()));

    // Populate cost and time matrices
    for (const auto& route : drivingInfo) {
        int cityA = cityIndex[route.city1];
        int cityB = cityIndex[route.city2];

        if (useHighway) {
            costMatrix[cityA][cityB] = route.highway_cost;
            timeMatrix[cityA][cityB] = route.highway_time;
        } else {
            costMatrix[cityA][cityB] = route.regular_cost;
            timeMatrix[cityA][cityB] = route.regular_time;
        }
    }

    // Initialize DP tables for cost and time
    vector<vector<double>> dpCost(1 << n, vector<double>(n, numeric_limits<double>::max()));
    vector<vector<double>> dpTime(1 << n, vector<double>(n, numeric_limits<double>::max()));
    vector<vector<int>> parentCost(1 << n, vector<int>(n, -1));
    vector<vector<int>> parentTime(1 << n, vector<int>(n, -1));

    dpCost[1][0] = 0;
    dpTime[1][0] = 0;

    // DP for cost and time
    for (int mask = 1; mask < (1 << n); ++mask) {
        for (int last = 0; last < n; ++last) {
            if (!(mask & (1 << last))) continue;

            for (int next = 0; next < n; ++next) {
                if (mask & (1 << next)) continue;

                int nextMask = mask | (1 << next);
                double newCost = dpCost[mask][last] + costMatrix[last][next];
                double newTime = dpTime[mask][last] + timeMatrix[last][next];

                // Update cost DP
                if (newCost < dpCost[nextMask][next]) {
                    dpCost[nextMask][next] = newCost;
                    parentCost[nextMask][next] = last;
                }

                // Update time DP
                if (newTime < dpTime[nextMask][next]) {
                    dpTime[nextMask][next] = newTime;
                    parentTime[nextMask][next] = last;
                }
            }
        }
    }

    // Retrieve the best path based on the chosen metric
    vector<CityVisit> visitList;
    double minValue = numeric_limits<double>::max();
    int lastCity = -1;

    // Determine whether to find min cost or min time
    for (int i = 1; i < n; ++i) {
        double valueWithReturn = (minimizeCost ? dpCost[allVisited][i] + costMatrix[i][0] : dpTime[allVisited][i] + timeMatrix[i][0]);
        
        if (valueWithReturn < minValue) {
            minValue = valueWithReturn;
            lastCity = i;
        }
    }

    // Backtrack to build the path
    int mask = allVisited;
    while (lastCity != -1) {
        // Retrieve the parent city index
        int parentCity = minimizeCost ? parentCost[mask][lastCity] : parentTime[mask][lastCity];

        // If there's a valid parent city, get the cost and time from drivingInfo
        double travelCost = 0.0;
        double travelTime = 0.0;

        // Find the travel cost and time from drivingInfo
        for (const auto& route : drivingInfo) {
            if (route.city1 == cities[parentCity] && route.city2 == cities[lastCity]) {
                travelCost = minimizeCost ? route.regular_cost : route.highway_cost;
                travelTime = minimizeCost ? route.regular_time : route.highway_time;
                break;  // Exit loop once we find the route
            }
        }

        // Create the visit record
        CityVisit visit = {
            cities[lastCity],
            travelCost,  // Cost to travel to current city
            travelTime,   // Time to travel to current city
            citycoordinates[cities[lastCity]]  // Get actual point from the citycoordinates map
        };

        visitList.push_back(visit);
        int temp = lastCity;
        lastCity = parentCity;  // Update lastCity to the parent
        mask ^= (1 << temp);  // Remove the visited city from the mask
    }

    reverse(visitList.begin(), visitList.end());
    return visitList;  // Return the selected path
}

double getTravelMetric(const vector<RouteInfo>& drivingInfo, const string& city1, const string& city2, bool useHighway, bool minimizeCost) {
    for (const auto& route : drivingInfo) {
        if ((route.city1 == city1 && route.city2 == city2) ||
            (route.city1 == city2 && route.city2 == city1)) {
            if (useHighway) {
                return minimizeCost ? route.highway_cost : route.highway_time;
            } else {
                return minimizeCost ? route.regular_cost : route.regular_time;
            }
        }
    }
    return INFINITY; // No path found between the cities
}

vector<CityVisit> tspDivideConquer(int left, int right, bool useHighway, bool minimizeCost, 
                                   vector<string>& cities, vector<RouteInfo>& drivingInfo, 
                                   unordered_map<string, cv::Point>& cityCoordinates) {
    vector<CityVisit> path;

    // Base case: only one city left (left == right)
    if (left == right) {
        string city = cities[left];
        cv::Point point = cityCoordinates[city];
        path.push_back({city, 0, 0, point});
        return path;
    }
    
    // Base case: only two cities left to connect
    if (right - left == 1) {
        double cost = 0, time = 0;
        string city1 = cities[left];
        string city2 = cities[right];
        cv::Point point1 = cityCoordinates[city1];
        cv::Point point2 = cityCoordinates[city2];

        // Look up cost and time in drivingInfo
        for (const auto& route : drivingInfo) {
            if (route.city1 == city1 && route.city2 == city2) {
                cost = useHighway ? route.highway_cost : route.regular_cost;
                time = useHighway ? route.highway_time : route.regular_time;
                break;
            }
        }

        // Push each city with its corresponding cost and time
        path.push_back({city1, 0, 0, point1});
        path.push_back({city2, cost, time, point2});
        return path;
    }

    // Recursive calls to divide cities into two groups
    int mid = left + (right - left) / 2;
    vector<CityVisit> leftPath = tspDivideConquer(left, mid, useHighway, minimizeCost, cities, drivingInfo, cityCoordinates);
    vector<CityVisit> rightPath = tspDivideConquer(mid + 1, right, useHighway, minimizeCost, cities, drivingInfo, cityCoordinates);

    // Ensure leftPath and rightPath are non-empty
    if (leftPath.empty() || rightPath.empty()) {
        cerr << "Error: Empty path encountered in recursion." << endl;
        return path;
    }

    // Merge paths by only adding the direct cost and time from the last city of leftPath to the first city of rightPath
    string lastCityLeft = leftPath.back().city;
    string firstCityRight = rightPath.front().city;
    double mergeCost = 0, mergeTime = 0;

    // Find route between last city of leftPath and first city of rightPath
    for (const auto& route : drivingInfo) {
        if ((route.city1 == lastCityLeft && route.city2 == firstCityRight) ) {
            mergeCost = useHighway ? route.highway_cost : route.regular_cost;
            mergeTime = useHighway ? route.highway_time : route.regular_time;
            break;
        }
    }

    // Only add the merge cost and time to the entry for the first city in the rightPath
    rightPath.front().cost = mergeCost;
    rightPath.front().time = mergeTime;

    // Merge the two halves into one final path
    leftPath.insert(leftPath.end(), rightPath.begin(), rightPath.end());
    return leftPath;

}

// Function to print the travel sequence
void printTravelSequence(const vector<CityVisit>& sequence) {
    double total_cost = 0, total_time = 0;
    for (const auto& step : sequence) {
        //cout << "City: " << step.city << ", Cost: " << step.cost << ", Time: " << step.time << step.point <<endl;
        total_cost += step.cost;
        total_time += step.time;
    }
    cout << "Total Cost: " << total_cost << ", Total Time: " << total_time << endl;
}

/// Function to write visitListR to CSV
void writeCSV(const string& filename, const vector<CityVisit>& visitList) {
    ofstream file(filename);
    if (file.is_open()) {
        file << "City,Cost,Time,PointX,PointY\n";  // CSV header
        for (const auto& visit : visitList) {
            file << visit.city << "," << visit.cost << "," << visit.time << ","
                 << visit.point.x << "," << visit.point.y << "\n";
        }
        file.close();
    } else {
        cerr << "Unable to open file: " << filename << endl;
    }
}

// Function to draw Initial Path without order (plot cities, draw paths, and annotate them with their index)
void drawInitialPath(cv::Mat& img, unordered_map<string, cv::Point> cityCoordinates)
{

    vector<cv::Point> city_list;
    for (const auto& city : cityCoordinates) {
        city_list.push_back(city.second); // Only push back the coordinates (cv::Point)
        cout << 7 << endl;
    }

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

void drawShortestGreedyPath(cv::Mat& img, const vector<CityVisit> visitListR)
{
    
    // Plot the cities and annotate them with their new greedy index
    for (size_t i = 0; i < visitListR.size(); i++)
    {
        // Plot city as a red circle
        cv::circle(img, visitListR[i].point, 10, cv::Scalar(0, 0, 255), -1);

        // Annotate the city with its new greedy index
        cv::putText(img, to_string(i), visitListR[i].point, cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw the greedy path with green lines based on the calculated path
    for (size_t i = 0; i < visitListR.size() - 1; i++)
    {
        cv::line(img, visitListR[i].point, visitListR[i+1].point, cv::Scalar(30, 150, 80), 4);
    }

    // Close the greedy path (draw line from last city to first city in path)
    if (!visitListR.empty())
    {
        cv::line(img, visitListR[visitListR.size() - 1].point, visitListR[0].point, cv::Scalar(30, 150, 80), 4);
    }
}

void drawShortestDpPath(cv::Mat& img, const vector<CityVisit> visitListR)
{
    
    // Plot the cities and annotate them with their new greedy index
    for (size_t i = 0; i < visitListR.size(); i++)
    {
        // Plot city as a red circle
        cv::circle(img, visitListR[i].point, 10, cv::Scalar(0, 0, 255), -1);

        // Annotate the city with its new greedy index
        cv::putText(img, to_string(i), visitListR[i].point, cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw the greedy path with green lines based on the calculated path
    for (size_t i = 0; i < visitListR.size() - 1; i++)
    {
        cv::line(img, visitListR[i].point, visitListR[i+1].point, cv::Scalar(170, 30, 220), 4);
    }

    // Close the greedy path (draw line from last city to first city in path)
    if (!visitListR.empty())
    {
        cv::line(img, visitListR[visitListR.size() - 1].point, visitListR[0].point, cv::Scalar(170, 30, 220), 4);
    }
}

void drawShortestDcPath(cv::Mat& img, const vector<CityVisit> visitListR)
{
    
    // Plot the cities and annotate them with their new greedy index
    for (size_t i = 0; i < visitListR.size(); i++)
    {
        // Plot city as a red circle
        cv::circle(img, visitListR[i].point, 10, cv::Scalar(0, 0, 255), -1);

        // Annotate the city with its new greedy index
        cv::putText(img, to_string(i), visitListR[i].point, cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0, 0, 0), 6);
    }

    // Draw the greedy path with green lines based on the calculated path
    for (size_t i = 0; i < visitListR.size() - 1; i++)
    {
        cv::line(img, visitListR[i].point, visitListR[i+1].point, cv::Scalar(20, 80, 50), 4);
    }

    // Close the greedy path (draw line from last city to first city in path)
    if (!visitListR.empty())
    {
        cv::line(img, visitListR[visitListR.size() - 1].point, visitListR[0].point, cv::Scalar(20, 80, 50), 4);
    }
}

/// Main processing function
void processPaths(vector<PathInfo>& paths, const vector<vector<CityVisit>>& visitLists, int algo_type) {
    for (size_t i = 0; i < paths.size(); ++i) {
        cout << paths[i].outputname << ": ";
        //print the total cost&time
        printTravelSequence(visitLists[i]);

        if (algo_type == 1){
            // Draw the path on the map
            drawShortestGreedyPath(paths[i].map, visitLists[i]);
        }
        else if(algo_type == 2){
            drawShortestDpPath(paths[i].map, visitLists[i]);
        }
        else if (algo_type == 3){
            drawShortestDcPath(paths[i].map, visitLists[i]);
        }
        
        // Save the image with the path drawn
        cv::imwrite(paths[i].path_name, paths[i].map);

        // Write the corresponding visitList to CSV
        writeCSV(paths[i].csv_filename, visitLists[i]);
    }
}

int main() {
    //type of algorithmn
    int type; // 1 for greedy, 2 for DP, 3 for DC.
    //file names for dataset
    string filenameD = "Dataset/Dataset_Driving1.csv";
    string filenameC = "Dataset/Cites.csv";
    string filenameCC = "Dataset/Dataset_Coordinate.csv";
    
    //load the data seperatly, initialize the vector that store the final path
    vector<string> cities = loadCities(filenameC);
    vector<RouteInfo> drivingInfo = loadRoutesInfo(filenameD);
    unordered_map<string, cv::Point> cityCoordinates = loadCityCoordinates(filenameCC);
    
    // Load the image (replace with your image path)
    string image_file = "Image/Europe.png";
    
    
    // Initialize the path info vector and sample image
    vector<PathInfo> greedy_paths = {
        {cv::imread("Image/Europe.png"),"Greedy Regular Road Min Cost Path" ,"./Image/E2Driving_Greedy_MinRegularCost_Path.png", "./Table/E2Greedy_Order_R_mincost.csv"},
        {cv::imread("Image/Europe.png"),"Greedy Regular Road Min Time Path" ,"./Image/E2Driving_Greedy_MinRegulartime_Path.png", "./Table/E2Greedy_Order_R_mintime.csv"},
        {cv::imread("Image/Europe.png"),"Greedy Highway Min Cost Path","./Image/E2Driving_Greedy_MinHighwaycost_Path.png", "./Table/E2Greedy_Order_H_mincost.csv"},
        {cv::imread("Image/Europe.png"),"Greedy Highway Min Time Path" ,"./Image/E2Driving_Greedy_MinHighwaytime_Path.png", "./Table/E2Greedy_Order_H_mintime.csv"},
    };

    //Greedy solution for both regular and highway
    vector<CityVisit> Greedy_visitList_R_Mincost = travelGreedy(cities, drivingInfo, cityCoordinates, false, true);
    vector<CityVisit> Greedy_visitList_R_Mintime = travelGreedy(cities, drivingInfo, cityCoordinates, false, false);
    vector<CityVisit> Greedy_visitList_H_Mincost = travelGreedy(cities, drivingInfo, cityCoordinates, true, true);
    vector<CityVisit> Greedy_visitList_H_Mintime = travelGreedy(cities, drivingInfo, cityCoordinates, true, false);

    vector<vector<CityVisit>> Greedy_visitLists = {
    Greedy_visitList_R_Mincost,
    Greedy_visitList_R_Mintime,
    Greedy_visitList_H_Mincost,
    Greedy_visitList_H_Mintime
    };

    // Process each path, saving images and writing CSV files
    type = 1;
    processPaths(greedy_paths, Greedy_visitLists, type);



    //Result for DP(time and cost for highway and regular road)
    // Initialize the path info vector and sample image
    vector<PathInfo> dp_paths = {
        {cv::imread("Image/Europe.png"),"DP Regular Road Min Cost Path" ,"./Image/E2Driving_DP_MinRegularCost_Path.png", "./Table/E2DP_Order_R_mincost.csv"},
        {cv::imread("Image/Europe.png"),"DP Regular Road Min Time Path" ,"./Image/E2Driving_DP_MinRegulartime_Path.png", "./Table/E2DP_Order_R_mintime.csv"},
        {cv::imread("Image/Europe.png"),"DP Highway Min Cost Path","./Image/E2Driving_DP_MinHighwaycost_Path.png", "./Table/E2DP_Order_H_mincost.csv"},
        {cv::imread("Image/Europe.png"),"DP Highway Min Time Path" ,"./Image/E2Driving_DP_MinHighwaytime_Path.png", "./Table/E2DP_Order_H_mintime.csv"},
    };

    //DP solution for both regular and highway
    vector<CityVisit> DP_visitList_R_Mincost = tspDP(cities, drivingInfo, cityCoordinates, false, true);
    vector<CityVisit> DP_visitList_R_Mintime = tspDP(cities, drivingInfo, cityCoordinates, false, false);
    vector<CityVisit> DP_visitList_H_Mincost = tspDP(cities, drivingInfo, cityCoordinates, true, true);
    vector<CityVisit> DP_visitList_H_Mintime = tspDP(cities, drivingInfo, cityCoordinates, true, false);

    vector<vector<CityVisit>> DP_visitLists = {
    DP_visitList_R_Mincost,
    DP_visitList_R_Mintime,
    DP_visitList_H_Mincost,
    DP_visitList_H_Mintime
    };

    // Process each path, saving images and writing CSV files
    type = 2;
    processPaths(dp_paths, DP_visitLists, type);



    //Result for DC(time and cost for highway and regular road)
    // Initialize the path info vector and sample image
    vector<PathInfo> dc_paths = {
        {cv::imread("Image/Europe.png"),"DC Regular Road Min Cost Path" ,"./Image/E2Driving_DC_MinRegularCost_Path.png", "./Table/E2DC_Order_R_mincost.csv"},
        {cv::imread("Image/Europe.png"),"DC Regular Road Min Time Path" ,"./Image/E2Driving_DC_MinRegulartime_Path.png", "./Table/E2DC_Order_R_mintime.csv"},
        {cv::imread("Image/Europe.png"),"DC Highway Min Cost Path","./Image/E2Driving_DC_MinHighwaycost_Path.png", "./Table/E2DC_Order_H_mincost.csv"},
        {cv::imread("Image/Europe.png"),"DC Highway Min Time Path" ,"./Image/E2Driving_DC_MinHighwaytime_Path.png", "./Table/E2DC_Order_H_mintime.csv"},
    };

    //DP solution for both regular and highway
    int left = 0;
    int right = cities.size() - 1;
    vector<CityVisit> DC_visitList_R_Mincost = tspDivideConquer(left, right, false, true, cities, drivingInfo, cityCoordinates);
    vector<CityVisit> DC_visitList_R_Mintime = tspDivideConquer(left, right, false, false, cities, drivingInfo, cityCoordinates);
    vector<CityVisit> DC_visitList_H_Mincost = tspDivideConquer(left, right, true, true, cities, drivingInfo, cityCoordinates);
    vector<CityVisit> DC_visitList_H_Mintime = tspDivideConquer(left, right, true, false, cities, drivingInfo, cityCoordinates);
    
    vector<vector<CityVisit>> DC_visitLists = {
        DC_visitList_R_Mincost,
        DC_visitList_R_Mintime,
        DC_visitList_H_Mincost,
        DC_visitList_H_Mintime
    }; 

    // Process each path, saving images and writing CSV files
    type = 3;
    processPaths(dc_paths, DC_visitLists, type);

    
   
    
    vector<PathInfo> all_paths = greedy_paths;
    all_paths.insert(all_paths.end(), dp_paths.begin(), dp_paths.end());
    all_paths.insert(all_paths.end(), dc_paths.begin(), dc_paths.end());

for (const auto& path_info : all_paths) {
#ifdef _WIN32
    system(("start " + path_info.path_name).c_str());
#elif __APPLE__
    system(("open " + path_info.path_name).c_str());
#elif __linux__
    system(("xdg-open " + path_info.path_name).c_str());
#endif
}
    return 0;
}