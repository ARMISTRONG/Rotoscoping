#ifndef QTANDOPENCV_H
#define QTANDOPENCV_H

#include <QtWidgets/QMainWindow>
#include "ui_qtandopencv.h"

#include <head.h>
#include <roto.h>


#include <QWidget>  
#include <QImage>  
#include <QTimer>     // ���òɼ����ݵļ��ʱ��
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

//���������Ƶ���л�����Ƶ
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
	string convertPath(string input);  //opencv�е�·��
	bool xyTest(int x, int y ,int width,int height);//���xy�Ƿ���������
	void curvePolyFit(vector<float> & x, vector<float> & y, int n, vector<float> & a, int m, vector<float> & dt);
	void tdma(vector<float> & x, vector<float> & a, vector<float> & b, vector<float> & c, vector<float> & d);
	float curvePolyFitValue(float t, vector<float> & a, int m, int n);//������Ϻ�����ֵ
	float bezierStepN(Point x0, Point x1, Point x2, Point x3);
	void drawBezier(Mat & frame, Point x0, Point x1, Point x2, Point x3,Scalar & color);//��resize֮���ͼƬ�ϻ�
	void drawCurve();//�����Ƶ�����֮�󻭵�resize֮���ͼƬ��
	void drawAllControlPointFlag();//���ɵ�ǰ�㵱ǰ֡
	void drawControlPointFlag(int CurrentControlPoint,int x, int y);//ֱ�ӻ�õ�����
	void drawAllControlPoint();
	void drawControlPoint(int x,int y);//�������ֱ�ӻ�ȡ������
	void updateAllControlPoint(int dx,int dy);//���µ�ǰ֡���п��Ƶ㣬����flag���������ת����ԭͼ�ߴ��ƽ����
	void updateControlPoint(int CurrentControlPoint, int dx,int dy);//���µ�ǰ֡��ǰ���Ƶ㣬����flag������������ź͹���֮�������
	void updateFirstControlPoint(int dx, int dy);//���µ�һ�������һ�����Ƶ㣨λ����ͬ��������������ź͹���֮�������
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
	void openVideo();      // ����Ƶ�ļ�
	void readResultText(); // ��ȡ���н���ļ�
	void readNextFrame();  // ��ȡ��һ֡��Ϣ  
	void readCurrentFrame();//��ȡ��ǰ֡
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
	Mat frameOriginal;	//��������ԭͼ
	Mat frame;			//resize֮���ͼƬ
	Mat ControlPointFlag;//����ʾ�����ڵ�ͼ�񲿷֣����һ���㲻�ǿ��Ƶ���ֵΪ-1
	Mat KeyFrameMark;//ֱ�۵���ʾ�ؼ�֡��λ��
	int FrameNum;//֡��
	int height ;//��ʾ�ĸ�
	int width ;	//��ʾ�Ŀ�
	int Height;	//ʵ�ʵĸ�
	int Width;	//ʵ�ʵĿ�
	int scaledHeight;//���ź�ĸ�
	int scaledWidth; //���ź�Ŀ�
	int hTranslation;//ˮƽ�϶�
	int vTranslation;//��ֱ�϶�
	double scale;//��ԭ�ߴ�����ű���
	int CurrentFrame;
	int CurrentLayer;
	int xPre, yPre;//�㱻�ƶ�֮ǰ��λ��
	int CurrentControlPoint;//��ǰҪ�޸�λ�õĿ��Ƶ����
	Interpolation InterpolationFlag;
	int PolynomialFitN;
	roto *pRoto;
protected:
	void mousePressEvent(QMouseEvent *e);       //--��갴���¼�  
	void mouseMoveEvent(QMouseEvent *e);    //--����ƶ��¼�  
	void mouseReleaseEvent(QMouseEvent *e); //--����ͷţ��ɿ����¼�  
	void mouseDoubleClickEvent(QMouseEvent *e); //--���˫���¼�
	void wheelEvent(QWheelEvent *e);//--�����¼�
	void keyPressEvent(QKeyEvent *e);//�����¼�
};

#endif // QTANDOPENCV_H
