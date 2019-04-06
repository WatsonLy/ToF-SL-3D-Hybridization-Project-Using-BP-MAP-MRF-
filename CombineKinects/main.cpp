#include <iostream>
#include <sstream>
#include <cstdio>
#include <algorithm>
#include <assert.h>
#include <cstring>
#include "opencv2/opencv.hpp"
#include <OpenNI.h>

#include "Restore.h"
#include "image.h"
#include "misc.h"
#include "pnmfile.h"
#include "Grabber.h"
#include "NtKinect.h"

void doJob() {
	NtKinect kinect;
	

	while (1) {
		
		kinect.setRGB();
		kinect.setDepth();
		kinect.setInfrared();


		//*Here starts James's attempts to create a depth frame 
		cv::Mat frame(kinect.depthImage);
		//frame.convertTo(frame, CV_8U, 255.0 / 65536.0);//255/65536 is the scaling factor for converting from 16bit to 8bit
		cv::imwrite("V2Frame.pgm", frame);
		//*///Here ends James's attempts to create a depth frame

		cv::imshow("rgb", kinect.rgbImage);
		cv::imshow("depth", kinect.depthImage);
		cv::imshow("infrared", kinect.infraredImage);

	
		/* -Blue stuff
		for (int y = 0; y < kinect.depthImage.rows; y++) {
			for (int x = 0; x < kinect.depthImage.cols; x++) {
				UINT16 d = kinect.depthImage.at<UINT16>(y, x);
				DepthSpacePoint dp; dp.X = x; dp.Y = y;
				ColorSpacePoint cp;
				HRESULT hr = kinect.coordinateMapper->MapDepthPointToColorSpace(dp, d, &cp);
				
				


				if (hr != S_OK) continue;
				if (d > 3000 || d < 500) {
					int ax = (int)cp.X;
					int ay = (int)cp.Y;
					cv::rectangle(kinect.rgbImage, cv::Rect(ax - 2, ay - 2, 4, 4), cv::Scalar(255, 0, 0), 2);
				}
			}
		}
		

		cv::imshow("rbg near", kinect.rgbImage);
		*/

		auto key = cv::waitKey(1);
		if (key == 'q') break;
	}
	cv::destroyAllWindows();
}

int main(int argc, char** argv) {

	try {
		
		//----------------------------------------------------------------------------------------------------------------
		//Grabbing Depth Frames
		//----------------------------------------------------------------------------------------------------------------
		//V1
		Grabber grabber;
		grabber.InitOpenNI();
		grabber.InitDevice();
		grabber.InitDepthStream();

		openni::VideoStream* streams[] = { grabber.depth_stream_ };
		int readyStream = -1;
		auto rc = openni::OpenNI::waitForAnyStream(streams, 1, &readyStream, 2000);
		if (rc != openni::STATUS_OK)
		{
			printf("Wait failed! (timeout is %d ms)\n%s\n", 2000, openni::OpenNI::getExtendedError());
		}
		cv::waitKey(50);
		cv::Mat DepthV1 = grabber.CapturePsenseDepthFrame();
		//cv::cvtColor(DepthV1, DepthV1, cv::COLOR_RGB2GRAY);

		//V2
		NtKinect kinectV2;

		for (int x = 0; x < 20; x++) {
			kinectV2.setDepth();
		}
		cv::Mat DepthV2(kinectV2.depthImage);
		
		//Shows the Depth Frames
		cv::imshow("KinectV1", DepthV1);
		cv::imshow("KinectV2", DepthV2);


		//----------------------------------------------------------------------------------------------------------------
		//Registering V2 to V1
		//----------------------------------------------------------------------------------------------------------------
		std::cout << "Registering Images from Kinects..." << endl;
		cv::Mat imReg, h, hDepth;
			
		hDepth = (cv::Mat_<float>(3, 3) << 1.5277089, -0.026776986, -0.0000358, -0.010839382,	1.5996544, -0.0001383637, -76.933472, -61.545399,	1);		
		
		hDepth = hDepth.t();

		warpPerspective(DepthV2, imReg, hDepth, DepthV1.size());
		cv::imshow("Registered KinectV2", imReg);


		//----------------------------------------------------------------------------------------------------------------
		//Cropping Images
		//----------------------------------------------------------------------------------------------------------------
		
		//int startX = 0, startY = 0, width = imReg.cols, height = imReg.rows;//Full
		//int startX = 0, startY = 0, width = imReg.cols, height = 360;		//Large
		int startX = 210, startY = 95, width = 225, height =265;			//Medium
		//int startX = 300, startY = 180, width = 55, height = 200;			//Small
		
		//int startX = 310, startY = 200, width = 35, height = 33;			//at depth 1.6
		//int startX = 311, startY = 206, width = 29, height = 29;			//at depth 1.9
		//int startX = 263, startY = 96, width = 120, height = 191;			//at depth 2.2
		//int startX = 217, startY = 135, width = 227, height = 203;		//Qualitative
		//int startX = 158, startY = 0, width = 290, height = 335;			//Multipath
		//int startX = 168, startY = 140, width = 100, height = 220;		//Multipath (Small)

		cv::Rect ROI(startX, startY, width, height);
		cv::Mat CropDepthV1 = DepthV1(ROI).clone();
		cv::Mat CropimReg = imReg(ROI).clone();
		
		double min, max;
		cv::minMaxIdx(CropDepthV1, &min, &max);
		cv::Mat ScaledV1, ScaledV2;
		cv::convertScaleAbs(CropDepthV1, ScaledV1, 255 / max);
		cv::convertScaleAbs(CropimReg, ScaledV2, 255 / max);

		cv::imshow("Cropped KinectV1", ScaledV1);
		cv::imshow("Cropped KinectV2", ScaledV2);

		//save images
		cv::imwrite("KinectV1.tif", DepthV1);
		cv::imwrite("KinectV2.tif", DepthV2);
		cv::imwrite("Registered KinectV2.tif", imReg);
		cv::imwrite("Cropped KinectV1.tif", CropDepthV1);
		cv::imwrite("Cropped KinectV2.tif", CropimReg);

		//----------------------------------------------------------------------------------------------------------------
		//Restoration using Belief Propogation
		//----------------------------------------------------------------------------------------------------------------
		cout << endl << "Press 'q' to start the restoration process using belief propagation" << endl;
		while (1) {
			auto key = cv::waitKey(1);
			if (key == 'q') break;
		}

		//Convert Mat object images to "image" objects.
		uint16_t* v1Arr = new uint16_t[CropDepthV1.rows*CropDepthV1.cols];
		uint16_t* v2Arr = new uint16_t[CropimReg.rows*CropimReg.cols];

		for (int i = 0; i < CropDepthV1.rows; ++i) {
			for (int j = 0; j < CropDepthV1.cols; ++j) {
				v1Arr[i*CropDepthV1.cols+j] = CropDepthV1.at<uint16_t>(i, j);
				//std::cout << (v1Arr[i*CropDepthV1.cols+j]);
				//std::cout << " ";
			}
		}
		
		for (int i = 0; i < CropimReg.rows; ++i) {
			for (int j = 0; j < CropimReg.cols; ++j) {
				v2Arr[i*CropimReg.cols+j] = CropimReg.at<uint16_t>(i, j);
				//std::cout << (v2Arr[i*CropimReg.cols+j]);
				//std::cout << " ";
			}
		}
		

		image<uint16_t> * inputV1 = new image<uint16_t>(CropDepthV1.cols, CropDepthV1.rows);
		inputV1->data = v1Arr;
		
		for (int i = 0; i < CropDepthV1.rows; ++i)
		{
			inputV1->access[i] = &inputV1->data[i*CropDepthV1.cols];
		}
			

		image<uint16_t> * inputV2 = new image<uint16_t>(CropimReg.cols, CropimReg.rows);
		inputV2->data = v2Arr;

		for (int i = 0; i < CropimReg.rows; ++i)
		{
			inputV2->access[i] = &inputV2->data[i*CropimReg.cols];
			
		}

		Restore SmartBot;
		image<uint16_t> *out;

		out = SmartBot.restore_ms(inputV1, inputV2);

		/*
		for (int i = 0; i < out->width() * out->height(); i++) {
			std::cout << out->data[i];
			std::cout << " ";
		}*/
		

		uint16_t* outArr = new uint16_t[CropDepthV1.rows*CropDepthV1.cols];
		outArr = out->data;

		cv::Mat restored = cv::Mat(CropDepthV1.rows, CropDepthV1.cols, CV_16UC1);
		for (int i = 0; i < CropDepthV1.rows; ++i) {
			for (int j = 0; j < CropDepthV1.cols; ++j) {
				restored.at<uint16_t>(i, j) = outArr[i*CropDepthV1.cols + j];
			}
		}

		cv::minMaxIdx(restored, &min, &max);
		cv::Mat adjMap;
		cv::convertScaleAbs(restored, adjMap, 255 / max);
		cv::imshow("Restored", adjMap);

		cv::imwrite("Restored.tif", restored);
		cv::imwrite("RestoredScaled.tif", adjMap);

		//delete v1Arr, v2Arr;
		//delete inputV1, inputV2;

		while (1) {
			auto key = cv::waitKey(1);
			if (key == 'q') break;
		}
	}
	catch (exception &ex) {
		cout << ex.what() << endl;
		string s;
		cin >> s;
	}
	return 0;
}