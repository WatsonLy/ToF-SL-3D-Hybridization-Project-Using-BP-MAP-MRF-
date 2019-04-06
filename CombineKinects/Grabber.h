#pragma once
#include <opencv2/opencv.hpp>
#include <OpenNI.h>

#include <iostream>

class Grabber
{
public:

	cv::Mat CapturePsenseDepthFrame();
	openni::VideoStream*   depth_stream_;

	void InitOpenNI();
	void InitDevice();
	void InitDepthStream();
	void InitIRStream();
	void InitColorStream();
	void Run();

	Grabber();
	~Grabber();

private:
	//cv::Mat CapturePsenseDepthFrame();
	void CapturePsenseDepthFrame2();
	cv::Mat CapturePsenseIRFrame();
	void CapturePsenseIRFrame2();
	void CapturePsenseColorFrame();
	cv::Mat ChangeDepthForDisplay(const cv::Mat&);

	openni::Device*        device_;
	//openni::VideoStream*   depth_stream_;
	openni::VideoStream*   ir_stream_;
	openni::VideoStream*   color_stream_;
	openni::VideoFrameRef* depth_frame_;
	openni::VideoFrameRef* ir_frame_;
	openni::VideoFrameRef* color_frame_;


};

