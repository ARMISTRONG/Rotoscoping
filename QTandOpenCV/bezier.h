#ifndef BEZIER_H
#define BEZIER_H
#pragma once
#include "head.h"


/*
void MouseCallbackFuncFirstFrame(int mouseevent,int x,int y,int flag, void * param);
void MouseCallbackFuncOtherFrame(int mouseevent,int x,int y,int flag, void * param);
void MouseCallbackFuncDisplacement(int mouseevent,int x,int y,int flag, void * param);



void MouseCallbackFuncFirstFrame(int mouseevent, int x, int y, int flag, void * param){

	Rotoscoping* RS = (Rotoscoping*)param;

	static bool Pressed = false;
	static bool IsControlPoint;
	static bool IsAdjustPoint;
	static int ControlPointId;//Id即为Flag
	static int AdjustPointId;
	static int xpre, ypre;
	static int CurPointId = 0;
	static int count = 0;//控制点计数，为4的倍数的时候画Bezier曲线

	if (mouseevent == CV_EVENT_LBUTTONDOWN)
	{
		Pressed = true;

		// If we press on a control point that already exists, then we move it
		ControlPointId = RS->GetControlPointFlag(x, y);
		if (ControlPointId >= 0){  // This a control point
			IsControlPoint = true;
			cout << "\t\t\t Updating Control Point " << ControlPointId << endl;

			//Erase the original Flag!!!! i.e, set to -1
			cvCircle(RS->ControlPointFlag, RS->ControlPoints[ControlPointId], FlagRadius, cvScalar(-1), -1);
		}
		else{
			IsControlPoint = false;
			++count;//控制点增加
		}

		cvCircle(RS->CurrentImageShow, cvPoint(x, y), ControlPointRadius, cvScalar(0, 0, 255), -1);

		xpre = x;
		ypre = y;

	}
	else if (mouseevent == CV_EVENT_LBUTTONUP)
	{
		Pressed = false;
		if (IsControlPoint){
			//Update the Control Point Flag in order to search!!!!
			cvCircle(RS->ControlPointFlag, cvPoint(x, y), FlagRadius, cvScalar(ControlPointId), -1);
		}
		else
		{
			//Bezier adjutment point (yellow)
			//cvCircle(RS->CurrentImageShow,cvPoint(x,y),ControlPointRadius,cvScalar(0,0,255),-1);   

			//Mark the Control Point in order to search!!!!
			RS->ControlPoints.push_back(cvPoint(xpre, ypre));
			cvCircle(RS->ControlPointFlag, cvPoint(xpre, ypre), FlagRadius, cvScalar(RS->ControlPoints.size() - 1), -1);  //  Mark the Flag in order to search!!!!

			RS->DrawCurve(1);
		}

	}
	else if (mouseevent == CV_EVENT_MOUSEMOVE)
	{
		if (Pressed)
		{
			if (IsControlPoint){
				RS->UpdateControlPoint(ControlPointId, x, y);
				cvCopy(RS->CurrentImage, RS->CurrentImageShow);
				//draw previous point
				RS->DrawCurve(1);
			}
			else{
				cvCopy(RS->CurrentImage, RS->CurrentImageShow);
				//draw previous point
				RS->DrawCurve(1);

				cvCircle(RS->CurrentImageShow, cvPoint(xpre, ypre), ControlPointRadius, cvScalar(0, 0, 255), -1);// control point that has not been push into the vector (red)


				int position = RS->ControlPoints.size();
				if (position >= 4 && position % 3 == 1){
					CvPoint BezierPoint[SAMPLESperCURVE];
					Bezier::GetBezierPointer()->GenerateCurve(RS->ControlPoints[position - 4], RS->ControlPoints[position - 3], RS->ControlPoints[position - 2], RS->ControlPoints[position - 1], SAMPLESperCURVE, BezierPoint);
					for (int i = 0; i<SAMPLESperCURVE; i++){
						cvCircle(RS->CurrentImageShow, BezierPoint[i], CurvePointRadius, cvScalar(0, 255, 255), -1); // Bezier points are yellow
					}
				}
			}

		}

	}
}


void MouseCallbackFuncOtherFrame(int mouseevent, int x, int y, int flag, void * param){
	Rotoscoping* RS = (Rotoscoping*)param;

	static bool Pressed = false;
	static bool IsControlPoint;
	static bool IsAdjustPoint;
	static int ControlPointId;
	static int AdjustPointId;
	static int xpre, ypre;
	static int CurPointId = 0;


	if (mouseevent == CV_EVENT_LBUTTONDOWN)
	{
		Pressed = true;

		// If we press on a control point that already exists, then we move it
		ControlPointId = RS->GetControlPointFlag(x, y);
		if (ControlPointId >= 0)  // This a control point
		{
			IsControlPoint = true;
			cout << "\t\t\t Updating Control Point " << ControlPointId << endl;
			//Erase the original Flag!!!!i.e,set to -1
			cvCircle(RS->ControlPointFlag, RS->ControlPoints[ControlPointId], FlagRadius, cvScalar(-1), -1);
			//新加的一部分
			if (ControlPointId == 0){
				cvCircle(RS->ControlPointFlag, RS->ControlPoints[ControlPointId + 1], FlagRadius, cvScalar(-1), -1);
			}
			else if (ControlPointId == RS->ControlPoints.size() - 1){
				cvCircle(RS->ControlPointFlag, RS->ControlPoints[ControlPointId - 1], FlagRadius, cvScalar(-1), -1);
			}
			else if (ControlPointId % 3 == 0){
				cvCircle(RS->ControlPointFlag, RS->ControlPoints[ControlPointId - 1], FlagRadius, cvScalar(-1), -1);
				cvCircle(RS->ControlPointFlag, RS->ControlPoints[ControlPointId + 1], FlagRadius, cvScalar(-1), -1);
			}
		}
		else{
			IsControlPoint = false;
		}
		xpre = x;
		ypre = y;

	}
	else if (mouseevent == CV_EVENT_LBUTTONUP)
	{
		Pressed = false;
		if (IsControlPoint){
			//Update the Control Point Flag in order to search!!!!
			cvCircle(RS->ControlPointFlag, cvPoint(x, y), FlagRadius, cvScalar(ControlPointId), -1);
			//新加的一部分
			double xTemp = RS->ControlPoints[ControlPointId].x;
			double yTemp = RS->ControlPoints[ControlPointId].y;
			if (ControlPointId == 0){
				cvCircle(RS->ControlPointFlag, cvPoint(RS->ControlPoints[ControlPointId + 1].x + x - xTemp, RS->ControlPoints[ControlPointId + 1].y + y - yTemp), FlagRadius, cvScalar(ControlPointId + 1), -1);
			}
			else if (ControlPointId == RS->ControlPoints.size() - 1){
				cvCircle(RS->ControlPointFlag, cvPoint(RS->ControlPoints[ControlPointId - 1].x + x - xTemp, RS->ControlPoints[ControlPointId - 1].y + y - yTemp), FlagRadius, cvScalar(ControlPointId - 1), -1);
			}
			else if (ControlPointId % 3 == 0){
				cvCircle(RS->ControlPointFlag, cvPoint(RS->ControlPoints[ControlPointId + 1].x + x - xTemp, RS->ControlPoints[ControlPointId + 1].y + y - yTemp), FlagRadius, cvScalar(ControlPointId + 1), -1);
				cvCircle(RS->ControlPointFlag, cvPoint(RS->ControlPoints[ControlPointId - 1].x + x - xTemp, RS->ControlPoints[ControlPointId - 1].y + y - yTemp), FlagRadius, cvScalar(ControlPointId - 1), -1);
			}
		}
	}
	else if (mouseevent == CV_EVENT_MOUSEMOVE)
	{
		if (Pressed)
		{
			if (IsControlPoint)
			{
				double xTemp = RS->ControlPoints[ControlPointId].x;
				double yTemp = RS->ControlPoints[ControlPointId].y;
				RS->UpdateControlPoint(ControlPointId, x, y);
				//新加一部分,使得控制点相邻的控制点一起移动
				if (ControlPointId == 0){
					RS->UpdateAdjacentControlPoint(ControlPointId + 1, x - xTemp, y - yTemp);
				}
				else if (ControlPointId == RS->ControlPoints.size() - 1){
					RS->UpdateAdjacentControlPoint(ControlPointId - 1, x - xTemp, y - yTemp);
				}
				else if (ControlPointId % 3 == 0){
					RS->UpdateAdjacentControlPoint(ControlPointId - 1, x - xTemp, y - yTemp);
					RS->UpdateAdjacentControlPoint(ControlPointId + 1, x - xTemp, y - yTemp);
				}

				cvCopy(RS->CurrentImage, RS->CurrentImageShow);
				//draw previous point
				RS->DrawCurve(1);
			}
		}

	}

}


void MouseCallbackFuncDisplacement(int mouseevent, int x, int y, int flag, void * param){
	Rotoscoping* RS = (Rotoscoping*)param;

	static bool Pressed = false;
	static bool IsControlPoint;
	static bool IsAdjustPoint;
	static int ControlPointId;
	static int AdjustPointId;
	static int xpre, ypre;
	static int CurPointId = 0;


	if (mouseevent == CV_EVENT_LBUTTONDOWN)
	{
		Pressed = true;
		RS->EraseAllPointFlag();
		xpre = x;
		ypre = y;
	}
	else if (mouseevent == CV_EVENT_LBUTTONUP)
	{
		Pressed = false;
		RS->UpdateAllPointFlag();

		cvCopy(RS->CurrentImage, RS->CurrentImageShow);//(src,dst,mask)
		//draw previous point
		RS->DrawCurve(1);
	}

	else if (mouseevent == CV_EVENT_MOUSEMOVE)
	{
		if (Pressed)
		{
			int dx = x - xpre;
			int dy = y - ypre;
			xpre = x;
			ypre = y;

			RS->DisplaceAllPoint(dx, dy);

			cvCopy(RS->CurrentImage, RS->CurrentImageShow);
			//draw previous point
			RS->DrawCurve(1);
		}

	}

}
*/



class Bezier{

Bezier(){}
Bezier(const Bezier& in){}
static Bezier * B;
public:
	

	static Bezier * GetBezierPointer(){
		return B;
	}

	void GenerateCurve(Point x0, Point x1, Point x2, Point x3,int IntermediatePointNum,vector<Point> & Output){
		//  Calculate  line1( determined by x0 and x1 )  and  line2 (determined by x2 and x3)
		double step = 1.0/(IntermediatePointNum+1);
		double t=step;	//  只求中间的点
		for (int i=0;i<IntermediatePointNum;i++,t+=step)
		{
			double coeff0 = (1-t)*(1-t)*(1-t);
			double coeff1 = 3*(1-t)*(1-t)*t;
			double coeff2 = 3*(1-t)*t*t;
			double coeff3 = t*t*t;
			
			Output[i].x = coeff0 * x0.x + coeff1 * x1.x + coeff2*x2.x + coeff3*x3.x;
			Output[i].y = coeff0 * x0.y + coeff1 * x1.y + coeff2*x2.y + coeff3*x3.y;
 		}
	}

};

Bezier* Bezier::B = new Bezier();



#endif