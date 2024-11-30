# Experiment 1

## Greedy Algorithm
- Start at the first city (or any chosen city).
- At each step, find the nearest unvisited city and travel to it.
- For each unvisited city, calculate the distance to the current city.
- Select the city with the minimum distance as the next city.
- Repeat the process** until all cities are visited.
- Return to the starting city to complete the loop.

> **Time Complexity**:  
> The greedy algorithm checks every unvisited city for the closest one, making it an **O(n²)** algorithm, where `n` is the number of cities.

> **Space Complexity**:  
> The space complexity is **O(n)**, where `n` is the number of cities, due to storing the list of visited cities and distances.

---

## Divide-and-Conquer Algorithm
- Divide: Split the list of cities into two halves based on their x-coordinates.
  - Sort the cities by x-coordinates to prepare for splitting.
  - Recursively divide the cities until each subset contains two or three cities.
- Conquer: Solve the problem for each half recursively.
  - For small subsets (with 2 or 3 cities), compute the shortest pair by brute-force.
  - Track the minimum distance found in both halves.
- Combine: After solving both halves, check for any shorter path between cities across the division line.
  - Specifically, check for the closest pair where one city is in the left half and the other is in the right half.

> **Time Complexity**:  
> The divide-and-conquer approach operates in **O(n log n)** time because it sorts the cities and divides them recursively.

> **Space Complexity**:  
> The space complexity is **O(n)** due to the recursive division and the space required to store city coordinates and distances.

---

## Dynamic Programming (TSP) Algorithm
- Use Dynamic Programming with Bitmasking to find the globally optimal solution for visiting all cities and returning to the starting point.
- Initialize a DP table `dp[mask][i]` where:
  - `mask` represents the set of visited cities as a bitmask.
  - `i` is the current city being visited.
  - `dp[mask][i]` stores the minimum distance to visit all cities in the set `mask` ending at city `i`.
- Set the base case: Starting from city 0, `dp[1][0] = 0` (starting at city 0 with only city 0 visited).
- For each subset of cities, calculate the minimum distance for visiting each city and update the DP table by considering all possible next steps:
  - Transition from city `u` to city `v` if city `v` has not yet been visited in the current subset.
  - Update `dp[mask | (1 << v)][v]` to keep track of the shortest path to visit `v`.
- Once all cities are visited, check the minimum cost to return to the starting city.
- Reconstruct the path by backtracking from the last visited city to the starting city.

> **Time Complexity**:  
> The dynamic programming algorithm for the TSP has a time complexity of **O(n² * 2ⁿ)** due to the bitmasking approach for subsets and the recursive path reconstruction.

> **Space Complexity**:  
> The space complexity is **O(n * 2ⁿ)** because we need to store the DP table that tracks all possible subsets of cities and their respective shortest paths.

---

## Comparison of Algorithms

- **Greedy Algorithm**: 
  - **Time Complexity**: **O(n²)**
  - **Space Complexity**: **O(n)**
  - Fast and easy to implement but may not produce the optimal solution.

- **Divide-and-Conquer**:
  - **Time Complexity**: **O(n log n)**
  - **Space Complexity**: **O(n)**
  - Effective for the closest-pair problem but is not typically used for solving the full TSP. It is more specialized for dividing the problem space and solving subproblems efficiently.

- **Dynamic Programming (TSP)**: 
  - **Time Complexity**: **O(n² * 2ⁿ)**
  - **Space Complexity**: **O(n * 2ⁿ)**
  - Guarantees the optimal solution but is computationally expensive for large datasets. This algorithm is the gold standard for solving the TSP exactly.


Number of Cities	Greedy (units)	Divide-and-Conquer (units)	Dynamic Programming (units)
18	10831.7	14487.3 units	8619.13
15	9752.78	10024.7	7987.9
12	8864.51	9752.9	7738.06
9	7275.01	9351.54	6488.12
6	6097.13	5839.92	5839.92

Number of Cities	Greedy (km)	Divide-and-Conquer (km)	Dynamic Programming (km)
18	9748.53 	13038.6 	7757.22 
15	8777.50 	9022.24 	7189.11 
12	7978.06 	8777.61 	6964.25 
9	6547.51 	8416.39 	5839.31 
6	5487.42 	5255.93 	5255.93 


# Experiment 2

## Greedy Algorithm
- Separate paths for minimizing cost and time for both regular roads and highways.
- At each step, find the next city that minimizes the respective metric (cost or time) among unvisited cities.
- Repeat the process until all cities are visited for each path type.

> **Time Complexity**> **:
Each iteration calculates the metric (cost or time) for every unvisited city, making it an O(n²) algorithm for each path type.

> **Space Complexity**:
The space complexity is O(n) due to storing visited cities and paths for each metric.

---

## Divide-and-Conquer Algorithm

- Divide: Split the cities into subsets recursively.
- Conquer: Solve for the shortest path (cost or time) within each subset.
- Combine: Find the best "bridge" connection between subsets to merge their paths, minimizing the chosen metric.

> **Time Complexity**:
O(n log n) due to recursive division and subset merging.

> **Space Complexity**:
O(n) to store intermediate paths and maintain visited status.

---

## Dynamic Programming (TSP) Algorithm
- Use DP with precomputed metrics (cost and time) for regular roads and highways.
- DP table dp[mask][i] tracks the minimum cost or time to visit all cities in mask ending at city i.
- Base case: Start with only the first city visited (dp[1][0] = 0).
- Iterative updates for each subset of cities:
- Transition between cities while minimizing the respective metric.
- Path reconstruction: Backtrack from the final state to construct the optimal path.

> **Time Complexity**:
O(n² * 2ⁿ) for each metric due to the subset enumeration and distance computation.

> **Space Complexity**:
O(n * 2ⁿ) for storing the DP table.


# Experiment 3

## Greedy Algorithm
- Separate paths for minimizing cost and time.
- At each step, find the next city that minimizes the respective metric (cost or time) among unvisited cities.
- Repeat the process until all cities are visited for each path type.

> **Time Complexity**:
Each iteration calculates the metric (cost or time) for every unvisited city, making it an O(n²) algorithm for each path type.

> **Space Complexity**:
The space complexity is O(n) due to storing visited cities and paths for each metric.

---

## Divide-and-Conquer Algorithm
- Divide: Split the cities into subsets recursively.
- Conquer: Solve for the shortest path (cost or time) within each subset.
- Combine: Find the best "bridge" connection between subsets to merge their paths, minimizing the chosen metric.

> **Time Complexity**:
O(n log n) due to recursive division and subset merging.

> **Space Complexity**:
O(n) to store intermediate paths and maintain visited status.

---

## Dynamic Programming (TSP) Algorithm
- Use DP with precomputed metrics (cost and time).
- DP table dp[mask][i] tracks the minimum cost or time to visit all cities in mask ending at city i.
- Base case: Start with only the first city visited (dp[1][0] = 0).
- Iterative updates for each subset of cities:
- Transition between cities while minimizing the respective metric.
- Path reconstruction: Backtrack from the final state to construct the optimal path.

> **Time Complexity**:
O(n² * 2ⁿ) for each metric due to the subset enumeration and distance computation.

> **Space Complexity**:
O(n * 2ⁿ) for storing the DP table.