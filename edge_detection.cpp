// Edge Detecion of Video using Canny Filter
// Author : Harmeet Singh
// Email : harmeet.esal@gmail.com
// Date : January 5, 2016

#include <opencv\cv.h> 
#include <opencv2/opencv.hpp>
#include <opencv\highgui.h>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;

int main() {
	Mat image;          // Generation of Matrix to store image
	Mat edges;			// Generation of Matrix to store the captured edges
	VideoCapture cap;          //Initiate Video capture
	cap.open(0);
	namedWindow("window", 1);          // Window to display input video
	namedWindow("fwindow", 1);		   // Window to display edges detected
	while (1) {
		cap >> image;          //Duplicate camera stream to "image"
		imshow("window", image);          //Display image
		waitKey(10);				// Delay of 10ms	
		Canny(image, edges, 10, 100);		// Canny filter applied to the image
		imshow("fwindow", edges);			// Display of detected edges
	}
	return 0;
}