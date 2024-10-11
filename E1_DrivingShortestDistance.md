# Algorithms for Finding the Shortest Path

## Greedy Algorithm for Shortest Path

- **Start** at the first city (or any chosen city).
- **At each step**, find the nearest unvisited city and travel to it.
- For each unvisited city, **calculate the distance** to the current city.
- **Select** the city with the minimum distance as the next city.
- **Repeat the process** until all cities are visited.
- **Return** to the starting city to complete the loop.

> **Time Complexity**:  
> The greedy algorithm checks every unvisited city for the closest one, making it an **O(n²)** algorithm, where `n` is the number of cities.

> **Space Complexity**:  
> The space complexity is **O(n)**, where `n` is the number of cities, due to storing the list of visited cities and distances.

> **Describe**:  
> This approach is intuitive and simple but does not guarantee the globally optimal solution. The resulting path may be suboptimal compared to more advanced algorithms.

---

## Divide-and-Conquer Algorithm for Shortest Path

- **Divide**: Split the list of cities into two halves based on their x-coordinates.
  - Sort the cities by x-coordinates to prepare for splitting.
  - Recursively divide the cities until each subset contains two or three cities.
- **Conquer**: Solve the problem for each half recursively.
  - For small subsets (with 2 or 3 cities), compute the shortest pair by brute-force.
  - Track the minimum distance found in both halves.
- **Combine**: After solving both halves, check for any shorter path between cities across the division line.
  - Specifically, check for the closest pair where one city is in the left half and the other is in the right half.

> **Time Complexity**:  
> The divide-and-conquer approach operates in **O(n log n)** time because it sorts the cities and divides them recursively.

> **Space Complexity**:  
> The space complexity is **O(n)** due to the recursive division and the space required to store city coordinates and distances.

> **Describe**:  
> The divide-and-conquer strategy helps reduce the complexity of finding the shortest path, especially when dealing with large datasets. However, it is mainly useful for finding the closest pairs of cities, rather than solving the full TSP problem.

---

## Dynamic Programming (TSP) Algorithm for Shortest Path

- **Use Dynamic Programming with Bitmasking** to find the globally optimal solution for visiting all cities and returning to the starting point.
- **Initialize a DP table** `dp[mask][i]` where:
  - `mask` represents the set of visited cities as a bitmask.
  - `i` is the current city being visited.
  - `dp[mask][i]` stores the minimum distance to visit all cities in the set `mask` ending at city `i`.
- **Set the base case**: Starting from city 0, `dp[1][0] = 0` (starting at city 0 with only city 0 visited).
- **For each subset of cities**, calculate the minimum distance for visiting each city and update the DP table by considering all possible next steps:
  - Transition from city `u` to city `v` if city `v` has not yet been visited in the current subset.
  - Update `dp[mask | (1 << v)][v]` to keep track of the shortest path to visit `v`.
- **Once all cities are visited**, check the minimum cost to return to the starting city.
- **Reconstruct the path** by backtracking from the last visited city to the starting city.

> **Time Complexity**:  
> The dynamic programming algorithm for the TSP has a time complexity of **O(n² * 2ⁿ)** due to the bitmasking approach for subsets and the recursive path reconstruction.

> **Space Complexity**:  
> The space complexity is **O(n * 2ⁿ)** because we need to store the DP table that tracks all possible subsets of cities and their respective shortest paths.

> **Describe**:  
> This algorithm guarantees the optimal solution to the TSP problem but has a higher computational complexity, making it impractical for large datasets.

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
18	10831.7	14260.4	8619.13
15	9752.78	10170.1	7987.9
12	8864.51	11857	7738.06
9	7275.01	11444.5	6488.12
6	6097.13	7861.85	5839.92

Number of Cities	Greedy (km)	Divide-and-Conquer (km)	Dynamic Programming (km)
18	9748.53 	12834.36 	7757.22 
15	8777.50 	9153.09 	7189.11 
12	7978.06 	10671.30 	6964.25 
9	6547.51 	10300.05 	5839.31 
6	5487.42 	7075.67 	5255.93 
Analysis of the Results:

## Vary Input Size(Number of Cities) for three algorithms
### Dynamic Programming Outperforms Other Algorithms:
As shown in Table 9, Dynamic Programming consistently provides the shortest distance across all input sizes. For example, with 18 cities, Dynamic Programming yields a total trip distance of 7757.22 km, which is significantly shorter than the Greedy (9748.53 km) and Divide-and-Conquer (12834.36 km) methods.
This trend holds for all other city counts (6 to 18), demonstrating the efficiency of Dynamic Programming in solving the Traveling Salesman Problem (TSP) optimally.
### Greedy Algorithm as a Moderate Approach:
The Greedy algorithm performs better than Divide-and-Conquer but is still outperformed by Dynamic Programming. For instance, with 15 cities, the Greedy algorithm results in a total trip distance of 8777.50 km, compared to 7189.11 km for Dynamic Programming.
Although Greedy's performance improves with fewer cities, it fails to find the shortest possible distance, as it prioritizes local optima rather than the global optimal solution.
### Divide-and-Conquer's Increasing Suboptimality:
Interestingly, the Divide-and-Conquer algorithm performs poorly compared to the other two, especially as the number of cities increases. For 18 cities, the Divide-and-Conquer method produces a total distance of 12834.36 km, far exceeding the other algorithms.
This suboptimality arises from its divide-and-conquer strategy, which, while efficient in some problem domains, struggles to maintain an optimal overall route when dealing with the global nature of the TSP.
Impact of Input Size:
For all algorithms, the total distance naturally increases as the number of cities grows. However, the rate of increase is different across algorithms. For example, the total distance for Greedy grows from 5487.42 km for 6 cities to 9748.53 km for 18 cities, while Divide-and-Conquer jumps from 7075.67 km to 12834.36 km for the same range. Dynamic Programming shows the smallest rate of increase, which reflects its efficiency in managing large input sizes.
