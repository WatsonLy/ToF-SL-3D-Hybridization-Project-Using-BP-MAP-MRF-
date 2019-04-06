#pragma once
#include <opencv2/opencv.hpp>
#include "opencv2/features2d.hpp"


class Register
{
public:

	void alignImages(cv::Mat &im1, cv::Mat &im2, cv::Mat &im1Reg, cv::Mat &h);

	Register();
	~Register();
};

