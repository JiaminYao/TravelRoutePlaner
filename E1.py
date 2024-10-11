# Generate Helicopter Shortest Distance Comparison of Different Algorithms
import matplotlib.pyplot as plt

# Extract the data from the table
cities = [18, 15, 12, 9, 6]
greedy = [9748.53, 8777.50, 7978.06, 6547.51, 5487.42]
divide_and_conquer = [12834.36, 9153.09, 10671.30, 10300.05, 7075.67]
dynamic_programming = [7757.22, 7189.11, 6964.25, 5839.31, 5255.93]
output_image_path = "./Image/Shortest_Comparison.png"

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(cities, greedy, marker='o', label='Greedy (km)', color='blue')
plt.plot(cities, divide_and_conquer, marker='o', label='Divide-and-Conquer (km)', color='red')
plt.plot(cities, dynamic_programming, marker='o', label='Dynamic Programming (km)', color='green')

plt.title('Helicopter Shortest Distance Comparison of Different Algorithms')
plt.xlabel('Number of Cities')
plt.ylabel('Total Distance (km)')
plt.grid(True)
plt.legend()

# Save the image
plt.savefig(output_image_path)
plt.close()

output_image_path