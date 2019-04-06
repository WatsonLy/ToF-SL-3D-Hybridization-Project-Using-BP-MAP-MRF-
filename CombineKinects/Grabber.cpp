///////////////////////////////////////////////////////////////////////////////
//
// Simple program that reads depth and color (RGB) images from Primensense
// camera using OpenNI2 and displays them using OpenCV.
//
// Ashwin Nanjappa
///////////////////////////////////////////////////////////////////////////////

#include "Grabber.h"
#include "NtKinect.h"

using namespace std;
using namespace cv;

void Grabber::InitOpenNI()
{
	auto rc = openni::OpenNI::initialize();
	if (rc != openni::STATUS_OK)
	{
		printf("Initialize failed\n%s\n", openni::OpenNI::getExtendedError());
		exit(0);
	}
}

void Grabber::InitDevice()
{
	device_ = new openni::Device();
	auto rc = device_->open(openni::ANY_DEVICE);
	if (rc != openni::STATUS_OK)
	{
		printf("Couldn't open device\n%s\n", openni::OpenNI::getExtendedError());
		exit(0);
	}
}

void Grabber::InitDepthStream()
{
	depth_stream_ = new openni::VideoStream();

	// Create depth stream from device
	if (device_->getSensorInfo(openni::SENSOR_DEPTH) != nullptr)
	{
		auto rc = depth_stream_->create(*device_, openni::SENSOR_DEPTH);
		if (rc != openni::STATUS_OK)
		{
			printf("Couldn't create depth stream\n%s\n", openni::OpenNI::getExtendedError());
			exit(0);
		}
	}

	// Get info about depth sensor
	const openni::SensorInfo& sensor_info = *device_->getSensorInfo(openni::SENSOR_DEPTH);
	const openni::Array<openni::VideoMode>& arr = sensor_info.getSupportedVideoModes();

	// Look for VGA mode in depth sensor and set it for depth stream
	for (int i = 0; i < arr.getSize(); ++i)
	{
		const openni::VideoMode& vmode = arr[i];
		if (vmode.getPixelFormat() == openni::PIXEL_FORMAT_DEPTH_1_MM &&
			vmode.getResolutionX() == 640 &&
			vmode.getResolutionY() == 480)
		{
			depth_stream_->setVideoMode(vmode);
			break;
		}
	}

	// Start the depth stream
	auto rc = depth_stream_->start();
	if (rc != openni::STATUS_OK)
	{
		printf("Couldn't start the depth stream\n%s\n", openni::OpenNI::getExtendedError());
		exit(0);
	}

	depth_frame_ = new openni::VideoFrameRef();
}

void Grabber::InitIRStream()
{
	ir_stream_ = new openni::VideoStream();

	// Create IR stream from device
	if (device_->getSensorInfo(openni::SENSOR_IR) != nullptr)
	{
		auto rc = ir_stream_->create(*device_, openni::SENSOR_IR);
		if (rc != openni::STATUS_OK)
		{
			printf("Couldn't create IR stream\n%s\n", openni::OpenNI::getExtendedError());
			exit(0);
		}
	}

	// Get info about IR sensor
	const openni::SensorInfo& sensor_info = *device_->getSensorInfo(openni::SENSOR_IR);
	const openni::Array<openni::VideoMode>& arr = sensor_info.getSupportedVideoModes();

	// Look for VGA mode in ir sensor and set it for ir stream
	for (int i = 0; i < arr.getSize(); ++i)
	{
		const openni::VideoMode& vmode = arr[i];
		if (vmode.getPixelFormat() == openni::PIXEL_FORMAT_RGB888 &&
			vmode.getResolutionX() == 640 &&
			vmode.getResolutionY() == 480)
		{
			ir_stream_->setVideoMode(vmode);
			break;
		}
	}

	// Start the IR stream
	auto rc = ir_stream_->start();
	if (rc != openni::STATUS_OK)
	{
		printf("Couldn't start the IR stream\n%s\n", openni::OpenNI::getExtendedError());
		exit(0);
	}

	ir_frame_ = new openni::VideoFrameRef();

}

void Grabber::InitColorStream()
{
	color_stream_ = new openni::VideoStream();

	if (device_->getSensorInfo(openni::SENSOR_COLOR) != nullptr)
	{
		auto rc = color_stream_->create(*device_, openni::SENSOR_COLOR);
		if (rc != openni::STATUS_OK)
		{
			printf("Couldn't create color stream\n%s\n", openni::OpenNI::getExtendedError());
			exit(0);
		}
	}

	// Get info about color sensor
	const openni::SensorInfo& sensor_info = *device_->getSensorInfo(openni::SENSOR_COLOR);
	const openni::Array<openni::VideoMode>& arr = sensor_info.getSupportedVideoModes();

	// Look for VGA mode and set it for color stream
	for (int i = 0; i < arr.getSize(); ++i)
	{
		const openni::VideoMode& vmode = arr[i];
		if (
			vmode.getResolutionX() == 1280 &&
			vmode.getResolutionY() == 960)
		{
			color_stream_->setVideoMode(vmode);
			break;
		}
	}

	// Note: Doing image registration earlier than this seems to fail
	if (device_->isImageRegistrationModeSupported(openni::IMAGE_REGISTRATION_DEPTH_TO_COLOR))
	{
		auto rc = device_->setImageRegistrationMode(openni::IMAGE_REGISTRATION_DEPTH_TO_COLOR);
		if (rc == openni::STATUS_OK)
			std::cout << "Depth to color image registration set success\n";
		else
			std::cout << "Depth to color image registration set failed\n";
	}
	else
	{
		std::cout << "Depth to color image registration is not supported!!!\n";
	}


	// Note: Doing image registration earlier than this seems to fail

	// Start color stream
	auto rc = color_stream_->start();
	if (rc != openni::STATUS_OK)
	{
		printf("Couldn't start the depth stream\n%s\n", openni::OpenNI::getExtendedError());
		exit(0);
	}

	color_frame_ = new openni::VideoFrameRef();
}

cv::Mat Grabber::CapturePsenseDepthFrame()
{
	auto rc = depth_stream_->readFrame(depth_frame_);
	if (rc != openni::STATUS_OK)
	{
		printf("Read failed!\n%s\n", openni::OpenNI::getExtendedError());
	}

	if (depth_frame_->getVideoMode().getPixelFormat() != openni::PIXEL_FORMAT_DEPTH_1_MM && depth_frame_->getVideoMode().getPixelFormat() != openni::PIXEL_FORMAT_DEPTH_100_UM)
	{
		printf("Unexpected frame format\n");
	}

	// Get pointer to Primesense depth frame
	openni::DepthPixel* dev_buf_ptr = (openni::DepthPixel*) depth_frame_->getData();

	// Copy frame data to OpenCV mat
	cv::Mat depth_mat(depth_frame_->getHeight(), depth_frame_->getWidth(), CV_16U, dev_buf_ptr);

	cv::Mat disp_mat = ChangeDepthForDisplay(depth_mat);

	//cv::imshow("Depth", disp_mat);
	return depth_mat;
}

void Grabber::CapturePsenseDepthFrame2()
{
	auto rc = depth_stream_->readFrame(depth_frame_);
	if (rc != openni::STATUS_OK)
	{
		printf("Read failed!\n%s\n", openni::OpenNI::getExtendedError());
	}

	if (depth_frame_->getVideoMode().getPixelFormat() != openni::PIXEL_FORMAT_DEPTH_1_MM && depth_frame_->getVideoMode().getPixelFormat() != openni::PIXEL_FORMAT_DEPTH_100_UM)
	{
		printf("Unexpected frame format\n");
	}

	// Get pointer to Primesense depth frame
	openni::DepthPixel* dev_buf_ptr = (openni::DepthPixel*) depth_frame_->getData();

	// Copy frame data to OpenCV mat
	cv::Mat depth_mat(depth_frame_->getHeight(), depth_frame_->getWidth(), CV_16U, dev_buf_ptr);

	cv::Mat disp_mat = ChangeDepthForDisplay(depth_mat);

	cv::imshow("Depth", disp_mat);
}

cv::Mat Grabber::CapturePsenseIRFrame()
{
	// Read from stream to frame
	auto rc = ir_stream_->readFrame(ir_frame_);
	if (rc != openni::STATUS_OK)
	{
		printf("Read failed!\n%s\n", openni::OpenNI::getExtendedError());
	}

	// Pointer to Primesense ir frame
	openni::PixelFormat* dev_buf_ptr = (openni::PixelFormat*) ir_frame_->getData();

	// Make mat from camera data
	cv::Mat ir_mat(ir_frame_->getHeight(), ir_frame_->getWidth(), CV_16U, dev_buf_ptr);
	// Convert to BGR format for OpenCV
	//cv::cvtColor(ir_mat, ir_mat, cv::COLOR_RGB2GRAY);

	//cv::imshow("IR", ir_mat);
	return ir_mat;
}

void Grabber::CapturePsenseIRFrame2()
{
	// Read from stream to frame
	auto rc = ir_stream_->readFrame(ir_frame_);
	if (rc != openni::STATUS_OK)
	{
		printf("Read failed!\n%s\n", openni::OpenNI::getExtendedError());
	}

	// Pointer to Primesense ir frame
	openni::RGB888Pixel* dev_buf_ptr = (openni::RGB888Pixel*) ir_frame_->getData();

	// Make mat from camera data
	cv::Mat ir_mat(ir_frame_->getHeight(), ir_frame_->getWidth(), CV_16U, dev_buf_ptr);
	// Convert to BGR format for OpenCV
	//cv::cvtColor(ir_mat, ir_mat, cv::COLOR_RGB2GRAY);

	cv::imshow("IR", ir_mat);

}

void Grabber::CapturePsenseColorFrame()
{
	// Read from stream to frame
	auto rc = color_stream_->readFrame(color_frame_);
	if (rc != openni::STATUS_OK)
	{
		printf("Read failed!\n%s\n", openni::OpenNI::getExtendedError());
	}

	// Pointer to Primesense color frame
	openni::RGB888Pixel* dev_buf_ptr = (openni::RGB888Pixel*) color_frame_->getData();

	// Make mat from camera data
	cv::Mat color_mat(color_frame_->getHeight(), color_frame_->getWidth(), CV_8UC3, dev_buf_ptr);
	// Convert to BGR format for OpenCV
	cv::cvtColor(color_mat, color_mat, cv::COLOR_RGB2GRAY);

	cv::imshow("Color", color_mat);
}

void Grabber::Run()
{
	//Kinect v1 init
	openni::VideoStream* streams[] = { depth_stream_, color_stream_ };
	cv::Mat DepthV1;

	//Kinect v2 init
	NtKinect kinectV2;

	int i = 0;
	while (1) {

		//V1
		int readyStream = -1;
		auto rc = openni::OpenNI::waitForAnyStream(streams, 2, &readyStream, 2000);
		if (rc != openni::STATUS_OK)
		{
			printf("Wait failed! (timeout is %d ms)\n%s\n", 2000, openni::OpenNI::getExtendedError());
		}
		cv::waitKey(50);
		DepthV1 = CapturePsenseDepthFrame();
		cv::imshow("KinectV1", DepthV1);

		//V2
		kinectV2.setDepth();
		cv::Mat DepthV2(kinectV2.depthImage);
		
		//resize v2 to be the size of v1
		cv::Mat frame2resized;
		cv::resize(DepthV2, frame2resized, cv::Size(DepthV1.size()), 0, 0, INTER_LINEAR);
		cv::imshow("KinectV2", frame2resized);

		auto key = cv::waitKey(1000);
		if (key == 'q') break;
	}
	
	//Saves frame from each
	cv::cvtColor(DepthV1, DepthV1, COLOR_RGB2GRAY);
	cv::imwrite("V1Frame.pgm", DepthV1);
	cv::Mat DepthV2(kinectV2.depthImage);
	cv::imwrite("V2Frame.pgm", DepthV2);
	
	/*
	int i = 0;
	while (1) {

		int readyStream = -1;
		auto rc = openni::OpenNI::waitForAnyStream(streams, 2, &readyStream, 2000);
		if (rc != openni::STATUS_OK)
		{
			printf("Wait failed! (timeout is %d ms)\n%s\n", 2000, openni::OpenNI::getExtendedError());
		}
		waitKey(1000);
		//for (;;) {
		kinectv1 = CapturePsenseDepthFrame();
		string name = "V1Frame" + to_string(i);
		string plate = ".pgm";
		string file = name + plate;
		cout << "\nConcat: " << file;

		printf("key pressed.\n");
		cv::imwrite(file, kinectv1);

		i++;


		
		int readyStream = -1;
		auto rc = openni::OpenNI::waitForAnyStream(streams, 2, &readyStream, 2000);
		if (rc != openni::STATUS_OK)
		{
			printf("Wait failed! (timeout is %d ms)\n%s\n", 2000, openni::OpenNI::getExtendedError());
			break;
		}
		//CapturePsenseIRFrame();
		//int i = 0;
		//string name = "Frame" + i;
		switch (readyStream)
		{
		case 0:
			CapturePsenseIRFrame2();
			//cv::imwrite(name + ".jpg", ir_frame_->_getFrame);
			break;
		case 1:
			CapturePsenseDepthFrame();
			break;
		default:
			printf("Unxpected stream\n");
		}

		char c = cv::waitKey(10);
		if ('q' == c)
			break;
			
	} */
	// Remove the IRstream and CapturePsenseIRFrame so Color can be run
	//delete(streams);
}

cv::Mat Grabber::ChangeDepthForDisplay(const cv::Mat& mat)
{
	assert(CV_16U == mat.type());

	const float depth_near = 500;
	const float depth_far = 5000;

	const float alpha = 255.0 / (depth_far - depth_near);
	const float beta = -depth_near * alpha;

	cv::Mat fmat;
	mat.convertTo(fmat, CV_32F);

	for (int r = 0; r < mat.rows; ++r)
	{
		for (int c = 0; c < mat.cols; ++c)
		{
			float v = fmat.at<float>(r, c) * alpha + beta;

			if (v > 255) v = 255;
			if (v < 0)   v = 0;

			fmat.at<float>(r, c) = v;
		}
	}

	cv::Mat bmat;
	fmat.convertTo(bmat, CV_8U);

	cv::Mat cmat;
	cv::cvtColor(bmat, cmat, cv::COLOR_GRAY2BGR);
	cv::applyColorMap(cmat, cmat, cv::COLORMAP_OCEAN);

	return cmat;
}

Grabber::Grabber()
{
}


Grabber::~Grabber()
{
}
