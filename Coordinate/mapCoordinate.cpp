#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>

// Global variables
cv::Mat img, img_display;
std::vector<cv::Point> click_positions;

// Mouse callback function
void onMouse(int event, int x, int y, int, void*)
{
    if (event == cv::EVENT_LBUTTONDOWN)
    {
        // Record the clicked position
        cv::Point clicked_point(x, y);
        click_positions.push_back(clicked_point);

        // Display the coordinates in the console
        std::cout << "Clicked at: (" << x << ", " << y << ")" << std::endl;

        // Annotate the image with a circle and the coordinates
        cv::circle(img_display, clicked_point, 10, cv::Scalar(0, 0, 255), -1); // Red circle
        cv::putText(img_display, std::to_string(x) + ", " + std::to_string(y),
                    cv::Point(x - 60, y - 30), cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(0, 0, 255), 2);

        // Display the updated image
        cv::imshow("Image", img_display);
    }
}

int main()
{
    // Path to the image
    std::string filename = "../Image/Europe.png";

    // Load the image
    img = cv::imread(filename);

    if (img.empty())
    {
        std::cerr << "Error: Unable to load image." << std::endl;
        return -1;
    }

    // Create a copy of the image for displaying annotations
    img_display = img.clone();

    // Display the image in a window
    cv::namedWindow("Image", cv::WINDOW_AUTOSIZE);
    cv::imshow("Image", img_display);

    // Set the mouse callback function
    cv::setMouseCallback("Image", onMouse);

    // Wait indefinitely for a key press
    cv::waitKey(0);

    return 0;
}
