#ifndef ROTO_H
#define ROTO_H
#include<head.h>



class Frame{
public:
	int FrameID;
	Mat  ControlPointFlag;  //If this pixel is not a Control Point, then its values is -1
	vector<Point> ControlPoints; // point specified by the user

	Frame(){}
	~Frame(){}
};

class layer{
public:
	vector<int> KeyframeIndex;
	vector<bool> KeyframeFlag;
	vector<Frame> frames;
	int FrameNum;
	int CurrentKeyframeIndex;
	bool FirstKeyframeDone;
	bool GetResultDone;

	layer(long iFrameNum):FrameNum(iFrameNum){
		CurrentKeyframeIndex = 0;
		vector<Frame> temp0(FrameNum);
		vector<bool> temp(FrameNum, false);
		frames = temp0;
		KeyframeFlag = temp;
		FirstKeyframeDone = false;
		GetResultDone = false;
	}

	~layer(){}

};

class roto{
public:
	long FrameNum;
	int Height;
	int Width;
	vector<layer> layers;

	roto(long iFrameNum,int iHeight,int iWidth):FrameNum(iFrameNum),Height(iHeight),Width(iWidth){
		layer a(iFrameNum);
		layers.push_back(a);
	}

	~roto(){}

};





#endif