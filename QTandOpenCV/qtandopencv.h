#ifndef QTANDOPENCV_H
#define QTANDOPENCV_H

#include <QtWidgets/QMainWindow>
#include "ui_qtandopencv.h"

#include <head.h>
#include <roto.h>


#include <QWidget>  
#include <QImage>  
#include <QTimer>     // 设置采集数据的间隔时间
#include <QFileInfo>
#include <QFileDialog>
#include <QMouseEvent>
#include <QTextCodec>
#include <QMessageBox>
#include <QProgressDialog>




#define FlagRadius  4
#define CurvePointRadius 1
#define ControlPointRadius 3
#define LineWidth 1
#define SAMPLESperCURVE 10
#define MaxPolyFitN 20
#define PI 3.1415926535





//interpolation
enum Interpolation
{
	Linear, 
	Parabolic,
	CubicSpline,
	Polynomial
};

//输入的是视频序列还是视频
enum InputType
{
	Video,
	Sequence
};




class QTandOpenCV : public QMainWindow
{
	Q_OBJECT

public:
	QTandOpenCV(QWidget *parent = 0);
	~QTandOpenCV();
	string convertPath(string input);  //opencv中的路径
	bool xyTest(int x, int y ,int width,int height);//检测xy是否在区域内
	void curvePolyFit(vector<float> & x, vector<float> & y, int n, vector<float> & a, int m, vector<float> & dt);
	void tdma(vector<float> & x, vector<float> & a, vector<float> & b, vector<float> & c, vector<float> & d);
	float curvePolyFitValue(float t, vector<float> & a, int m, int n);//计算拟合函数的值
	float bezierStepN(Point x0, Point x1, Point x2, Point x3);
	void drawBezier(Mat & frame, Point x0, Point x1, Point x2, Point x3,Scalar & color);//在resize之后的图片上画
	void drawCurve();//将控制点缩放之后画到resize之后的图片上
	void drawAllControlPointFlag();//生成当前层当前帧
	void drawControlPointFlag(int CurrentControlPoint,int x, int y);//直接获得的坐标
	void drawAllControlPoint();
	void drawControlPoint(int x,int y);//输入的是直接获取的坐标
	void updateAllControlPoint(int dx,int dy);//更新当前帧所有控制点，包括flag，输入的是转换到原图尺寸的平移量
	void updateControlPoint(int CurrentControlPoint, int dx,int dy);//更新当前帧当前控制点，包括flag，输入的是缩放和滚动之后的坐标
	void updateFirstControlPoint(int dx, int dy);//更新第一个和最后一个控制点（位置相同），输入的是缩放和滚动之后的坐标
	void updateResultLinear();
	void updateResult();
	void copyControlPoint();
	void showFrame(Mat input);
	void showKeyFrameMark();
	void drawKeyFrameMark();
	void createFile(string Path);
	void deleteFile(string Path);
	void saveResultText(string Path);
	void drawBezierInMask(Mat & frame, Point x0, Point x1, Point x2, Point x3,int type);
	Mat getMask(int iCurrentLayer, int iCurrentFrame);
	Mat getMaskBGRA(int iCurrentLayer, int iCurrentFrame);
	Point getCenter(vector<Point> ControlPoints);
	void scaleControlPoints(vector<Point> & ControlPoints, float scale);
	void rotateControlPoints(vector<Point> & ControlPoints,float angle);
	long long imageTypeInSequence(string path, struct _finddata_t &fileinfo);
	void getFiles(string path, vector<string>& files);
	string sequenceToVideo(string path);
private slots:
	void openVideo();      // 打开视频文件
	void readResultText(); // 读取已有结果文件
	void readNextFrame();  // 读取下一帧信息  
	void readCurrentFrame();//读取当前帧
	void addKeyframe();
	void addLayer();
	void deleteLayer();
	void setLayer();
	void lastKeyframe();
	void nextKeyframe();
	void saveResult();
	void getResult();
	void setLinearFit();
	void setParabolicFit();
	void setCubicSplineFit();
	void setPolynomialFit();
	void setPolynomialFitN();
	void hMove(int);
	void vMove(int);
	void zoomInScale();
	void zoomOutScale();
private:
	Ui::QTandOpenCVClass ui;
	InputType inputType;
	string videoPath;
	string videoName;
	VideoCapture video;
	vector<string> files;
	Mat frameOriginal;	//读进来的原图
	Mat frame;			//resize之后的图片
	Mat ControlPointFlag;//在显示区域内的图像部分，如果一个点不是控制点则值为-1
	Mat KeyFrameMark;//直观地显示关键帧的位置
	int FrameNum;//帧数
	int height ;//显示的高
	int width ;	//显示的宽
	int Height;	//实际的高
	int Width;	//实际的宽
	int scaledHeight;//缩放后的高
	int scaledWidth; //缩放后的宽
	int hTranslation;//水平拖动
	int vTranslation;//垂直拖动
	double scale;//和原尺寸的缩放比例
	int CurrentFrame;
	int CurrentLayer;
	int xPre, yPre;//点被移动之前的位置
	int CurrentControlPoint;//当前要修改位置的控制点序号
	Interpolation InterpolationFlag;
	int PolynomialFitN;
	roto *pRoto;
protected:
	void mousePressEvent(QMouseEvent *e);       //--鼠标按下事件  
	void mouseMoveEvent(QMouseEvent *e);    //--鼠标移动事件  
	void mouseReleaseEvent(QMouseEvent *e); //--鼠标释放（松开）事件  
	void mouseDoubleClickEvent(QMouseEvent *e); //--鼠标双击事件
	void wheelEvent(QWheelEvent *e);//--滚轮事件
	void keyPressEvent(QKeyEvent *e);//键盘事件
};

#endif // QTANDOPENCV_H
