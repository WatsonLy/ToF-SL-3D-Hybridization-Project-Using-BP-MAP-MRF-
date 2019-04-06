#include "Register.h"
#include "opencv2/opencv.hpp"
#include "opencv2/features2d.hpp"

const int MAX_FEATURES = 500;
const float GOOD_MATCH_PERCENT = 0.15f;

using namespace cv;

void Register::alignImages(cv::Mat &im1, cv::Mat &im2, cv::Mat &im1Reg, cv::Mat &h)
{
	
	//convert mat type
	//std::cout << im1.type();
	//std::cout << im2.type();
	//im1.convertTo(im1, CV_32F, 65536 / 255.0);//65536/255 is the scaling factor for converting from 8bit to 16bit
	//imshow("im1", im1);

	//change size of image 2
	//Mat im2Resized;
	//cv::resize(im2, im2Resized, cv::Size(im1.size()), 0, 0, INTER_LINEAR);

	// Convert images to grayscale
	Mat im1Gray, im2Gray;
	cvtColor(im1, im1Gray, COLOR_BGR2GRAY);
	cvtColor(im2, im2Gray, COLOR_BGR2GRAY);

	// Variables to store keypoints and descriptors
	std::vector<KeyPoint> keypoints1, keypoints2;
	Mat descriptors1, descriptors2;

	// Detect ORB features and compute descriptors.
	Ptr<Feature2D> orb = ORB::create(MAX_FEATURES);
	orb->detectAndCompute(im1Gray, Mat(), keypoints1, descriptors1);
	orb->detectAndCompute(im2Gray, Mat(), keypoints2, descriptors2);

	// Match features.
	std::vector<DMatch> matches;
	Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce"); // BruteForce-Hamming
	matcher->match(descriptors1, descriptors2, matches, Mat());

	// Sort matches by score
	std::sort(matches.begin(), matches.end());

	// Remove not so good matches
	const int numGoodMatches = matches.size() * GOOD_MATCH_PERCENT;
	matches.erase(matches.begin() + numGoodMatches, matches.end());


	// Draw top matches
	Mat imMatches;
	drawMatches(im1, keypoints1, im2, keypoints2, matches, imMatches);
	imwrite("matches.jpg", imMatches);


	// Extract location of good matches
	std::vector<Point2f> points1, points2;

	for (size_t i = 0; i < matches.size(); i++)
	{
		points1.push_back(keypoints1[matches[i].queryIdx].pt);
		points2.push_back(keypoints2[matches[i].trainIdx].pt);
	}

	
	// Find homography
	h = findHomography(points1, points2, RANSAC); // LMEDS, RANSAC, RHOS
	std::cout << "\n\n" << h << std::endl;


	// Use homography to warp image
	//warpPerspective(im1, im1Reg, h, im1.size());
	warpPerspective(im1, im1Reg, h, im1.size());

}


Register::Register()
{
}


Register::~Register()
{
}
