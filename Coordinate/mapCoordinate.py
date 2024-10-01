import numpy as np
import cv2
import matplotlib.pyplot as plt

def onclick(event):
    if event.xdata is not None and event.ydata is not None:
        x, y = int(event.xdata), int(event.ydata)
        print(f'Clicked at: ({x}, {y})')
        click_positions.append((x, y))
        
        text_x = x - 60
        text_y = y - 30
        # Display the coordinates on the image
        cv2.circle(img_rgb, (x, y), radius=10, color=(0, 0, 255), thickness=-1)
        cv2.putText(img_rgb, f'{x}, {y}', (text_x, text_y), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)
        plt.imshow(img_rgb)
        plt.axis('on')
        plt.show()

if __name__ == "__main__":
    # Path to the image
    filename = "/Users/yaojiamin/Documents/SIUE/Course/CS456 Algorithms/project/Europe.png"

    img = cv2.imread(filename)

    if img is None:
        print("Error: Unable to load image.")
    else:
        # Convert the image from BGR (OpenCV format) to RGB (Matplotlib format)
        img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

        # Global variable to store click positions
        click_positions = []

        # Display the image with Matplotlib and connect the click event
        fig, ax = plt.subplots()
        ax.imshow(img_rgb)
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.axis('on')
        plt.show()
