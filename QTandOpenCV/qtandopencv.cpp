#include "qtandopencv.h"

QTandOpenCV::QTandOpenCV(QWidget *parent)
: QMainWindow(parent)
{
	ui.setupUi(this);
	pRoto = NULL; 
	height = ui.label->height();
	width = ui.label->width();
	hTranslation = 0;
	vTranslation = 0;
	xPre = -1;
	yPre = -1;
	ControlPointFlag = Mat(height, width, CV_32SC1, Scalar(-1));//显示窗口的大小
	CurrentControlPoint = -1;
	InterpolationFlag = Linear;
	PolynomialFitN = 3;
	ui.npolynomialfit->setRange(3, 3);
	

	/*信号和槽*/
	connect(ui.open, SIGNAL(clicked()), this, SLOT(openVideo()));
	connect(ui.loadresult, SIGNAL(clicked()), this, SLOT(readResultText()));
	connect(ui.horizontalSlider, SIGNAL(valueChanged(int)), this, SLOT(readCurrentFrame()));
	connect(ui.keyframe, SIGNAL(clicked()), this, SLOT(addKeyframe()));
	connect(ui.addlayer, SIGNAL(clicked()), this, SLOT(addLayer()));
	connect(ui.deletelayer, SIGNAL(clicked()), this, SLOT(deleteLayer())); 
	connect(ui.setlayer, SIGNAL(valueChanged(int)), this, SLOT(setLayer()));
	connect(ui.lastkeyframe, SIGNAL(clicked()), this, SLOT(lastKeyframe()));
	connect(ui.nextkeyframe, SIGNAL(clicked()), this, SLOT(nextKeyframe()));
	connect(ui.saveresult, SIGNAL(clicked()), this, SLOT(saveResult()));
	connect(ui.getresult, SIGNAL(clicked()), this, SLOT(getResult()));
	connect(ui.linearfit, SIGNAL(clicked()), this, SLOT(setLinearFit()));
	connect(ui.parabolicfit, SIGNAL(clicked()), this, SLOT(setParabolicFit()));
	connect(ui.cubicsplinefit, SIGNAL(clicked()), this, SLOT(setCubicSplineFit()));
	connect(ui.polynomialfit, SIGNAL(clicked()), this, SLOT(setPolynomialFit()));
	connect(ui.npolynomialfit, SIGNAL(valueChanged(int)), this, SLOT(setPolynomialFitN()));
	connect(ui.horizontalBar, SIGNAL(valueChanged(int)), this, SLOT(hMove(int)));
	connect(ui.verticalBar, SIGNAL(valueChanged(int)), this, SLOT(vMove(int)));
	connect(ui.zoomIn, SIGNAL(clicked()), this, SLOT(zoomInScale()));
	connect(ui.zoomOut, SIGNAL(clicked()), this, SLOT(zoomOutScale()));
}

void QTandOpenCV::showKeyFrameMark(){
	drawKeyFrameMark();
	//QImage(uchar *data, int width, int height, int bytesPerLine, Format format, QImageCleanupFunction cleanupFunction = 0, void *cleanupInfo = 0);
	QImage image = QImage(KeyFrameMark.data, KeyFrameMark.cols, KeyFrameMark.rows, KeyFrameMark.step, QImage::Format_RGB888).rgbSwapped();//指定每行的字节数，不然会出错
	ui.keyframemark->setPixmap(QPixmap::fromImage(image));  // 将图片显示到label上  
}

void QTandOpenCV::drawKeyFrameMark(){
	KeyFrameMark = Mat(20, 1296, CV_8UC3, Scalar(100, 100, 100));//1296*20
	int KeyFrameSize = pRoto->layers[CurrentLayer].KeyframeIndex.size();
	double step = 1286.0 / (pRoto->FrameNum - 1);
	for (int i = 0; i < KeyFrameSize;++i){
		int x = step*pRoto->layers[CurrentLayer].KeyframeIndex[i] + 4 + 0.5;
		line(KeyFrameMark, Point(x, 0), Point(x, 16), Scalar(0, 0, 0), 1);
	}
	//QImage(uchar *data, int width, int height, int bytesPerLine, Format format, QImageCleanupFunction cleanupFunction = 0, void *cleanupInfo = 0);
	QImage image = QImage(KeyFrameMark.data, KeyFrameMark.cols, KeyFrameMark.rows, KeyFrameMark.step, QImage::Format_RGB888).rgbSwapped();//指定每行的字节数，不然会出错
	ui.keyframemark->setPixmap(QPixmap::fromImage(image));  // 将图片显示到label上  
}

bool QTandOpenCV::xyTest(int x, int y,int width,int height){
	if (x >= 0 && x < width &&y >= 0 && y < height) return true;
	else return false;
}

string QTandOpenCV::convertPath(string input){
	for (int i = 0; i < input.size();++i){
		if (input[i] == '\\') input[i] = '/';
	}
	return input;
}

//解三对角矩阵方程
void QTandOpenCV::tdma(vector<float> & x, vector<float> & a, vector<float> & b, vector<float> & c, vector<float> & d)//解三对角矩阵方程
{
	int N = x.size();
	int n;

	c[0] = c[0] / b[0];
	d[0] = d[0] / b[0];

	for (n = 1; n < N; n++) {
		float m = 1.0f / (b[n] - a[n] * c[n - 1]);
		c[n] = c[n] * m;
		d[n] = (d[n] - a[n] * d[n - 1]) * m;
	}
	x[N - 1] = d[N - 1];
	for (n = N - 2; n >= 0; n--){
		x[n] = d[n] - c[n] * x[n + 1];
	}
}


/******************************************
//参考 《常用算法程序集 （C语言描述 第三版)》
//最小二乘法
//x[n] y[n] 已知输入
//n输入点个数
//a[m] 返回m-1次拟合多项式的m个系数
//m  拟合多项式的项数,即拟合多项式的最高次为m-1
//dt[3]:
//dt[0]返回拟合多项式与各数据点误差的平方和
//dt[1]返回拟合多项式与各数据点误差的绝对值之和
//dt[2]返回拟合多项式与各数据点误差的绝对值的最大值
//拟合多项式的输出
//Y(x) = a0 + a1(x-X) + a2(x-X)^2 + …… am(x-X)^m
//其中X为已知点x的平均值
******************************************/
void QTandOpenCV::curvePolyFit(vector<float> & x, vector<float> & y, int n, vector<float> & a, int m, vector<float> & dt) //x[],y[],a[],dt[];
{
	int i, j, k;
	float z, p, c, g, q, d1, d2, s[20], t[20], b[20];
	for (i = 0; i <= m - 1; i++) a[i] = 0.0;
	if (m>n) m = n;
	if (m>20) m = 20;
	z = 0.0;
	for (i = 0; i <= n - 1; i++) z = z + x[i] / (1.0*n);
	b[0] = 1.0; d1 = 1.0*n; p = 0.0; c = 0.0;
	for (i = 0; i <= n - 1; i++)
	{
		p = p + (x[i] - z); c = c + y[i];
	}
	c = c / d1; p = p / d1;
	a[0] = c*b[0];
	if (m>1)
	{
		t[1] = 1.0; t[0] = -p;
		d2 = 0.0; c = 0.0; g = 0.0;
		for (i = 0; i <= n - 1; i++)
		{
			q = x[i] - z - p; d2 = d2 + q*q;
			c = c + y[i] * q;
			g = g + (x[i] - z)*q*q;
		}
		c = c / d2; p = g / d2; q = d2 / d1;
		d1 = d2;
		a[1] = c*t[1]; a[0] = c*t[0] + a[0];
	}
	for (j = 2; j <= m - 1; j++)
	{
		s[j] = t[j - 1];
		s[j - 1] = -p*t[j - 1] + t[j - 2];
		if (j >= 3)
		for (k = j - 2; k >= 1; k--)
			s[k] = -p*t[k] + t[k - 1] - q*b[k];
		s[0] = -p*t[0] - q*b[0];
		d2 = 0.0; c = 0.0; g = 0.0;
		for (i = 0; i <= n - 1; i++)
		{
			q = s[j];
			for (k = j - 1; k >= 0; k--)
				q = q*(x[i] - z) + s[k];
			d2 = d2 + q*q; c = c + y[i] * q;
			g = g + (x[i] - z)*q*q;
		}
		c = c / d2; p = g / d2; q = d2 / d1;
		d1 = d2;
		a[j] = c*s[j]; t[j] = s[j];
		for (k = j - 1; k >= 0; k--)
		{
			a[k] = c*s[k] + a[k];
			b[k] = t[k]; t[k] = s[k];
		}
	}
	dt[0] = 0.0; dt[1] = 0.0; dt[2] = 0.0;
	for (i = 0; i <= n - 1; i++)
	{
		q = a[m - 1];
		for (k = m - 2; k >= 0; k--)
			q = a[k] + q*(x[i] - z);
		p = q - y[i];
		if (fabs(p)>dt[2]) dt[2] = fabs(p);
		dt[0] = dt[0] + p*p;
		dt[1] = dt[1] + fabs(p);
	}
	return;
}

float QTandOpenCV::curvePolyFitValue(float t, vector<float> & a, int m, int n){
	float result = 0.0;
	result += a[0];
	float T = (n - 1) / 2.0;//t的均值，0~n-1
	for (int i = 1; i<m; ++i){
		if ((t - T) != 0) result += a[i] * pow((t - T), i);
	}
	return result;
}


float QTandOpenCV::bezierStepN(Point x0, Point x1, Point x2, Point x3){
	float Vx = x2.x - 2 * x1.x + x0.x;
	float Wx = x3.x - 3 * (x2.x - x1.x) - x0.x;
	float Tx = -Vx / Wx;
	float Nx = max(abs(x3.x - x2.x), abs(x1.x - x0.x));
	if (Tx >= 0 && Tx <= 1){
		Nx = 3 * max(Nx, fabs(x1.x - x0.x + Tx*Vx));
	}
	else{
		Nx *= 3;
	}

	float Vy = x2.y - 2 * x1.y + x0.y;
	float Wy = x3.y - 3 * (x2.y - x1.y) - x0.y;
	float Ty = -Vy / Wy;
	float Ny = max(abs(x3.y - x2.y), abs(x1.y - x0.y));
	if (Ty >= 0 && Ty <= 1){
		Ny = 3 * max(Ny, fabs(x1.y - x0.y + Ty*Vy));
	}
	else{
		Ny *= 3;
	}

	return max(Nx, Ny);
}

//没有反走样的bezier曲线绘制
void QTandOpenCV::drawBezier(Mat & frame, Point x0, Point x1, Point x2, Point x3,Scalar & color){
	int BezierPointSize = bezierStepN(x0, x1, x2, x3);
	double step = 1.0 / (BezierPointSize + 1);

	double t = step;	//  只求中间的点
	double xLast = -1,yLast=-1;
	for (int i = 0; i<BezierPointSize; i++, t += step)
	{
		double coeff0 = (1 - t)*(1 - t)*(1 - t);
		double coeff1 = 3 * (1 - t)*(1 - t)*t;
		double coeff2 = 3 * (1 - t)*t*t;
		double coeff3 = t*t*t;
		double x = coeff0 * x0.x + coeff1 * x1.x + coeff2*x2.x + coeff3*x3.x;
		double y = coeff0 * x0.y + coeff1 * x1.y + coeff2*x2.y + coeff3*x3.y;

		//取最近点
		int xTemp = x + 0.5;
		int yTemp = y + 0.5;
		if (xTemp == xLast && yTemp == yLast){
			continue;
		}
		if (xyTest(xTemp, yTemp, frame.cols, frame.rows)){
			for (int j = 0; j < 3; ++j){
				frame.at<Vec3b>(yTemp, xTemp)[j] = color[j];
			}
		}
		xLast = xTemp;
		yLast = yTemp;
		
	}
}

void QTandOpenCV::drawCurve(){
	vector<Point> ControlPointTemp(pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints);
	for (int i = 0; i < ControlPointTemp.size(); ++i){//缩放到显示大小,并对齐显示区域
		ControlPointTemp[i].x = ControlPointTemp[i].x *scale;
		ControlPointTemp[i].y = ControlPointTemp[i].y *scale;
	}
	if (ControlPointTemp.size()>1) 
		line(frame, ControlPointTemp[0], ControlPointTemp[1], Scalar(255, 0, 0), 1);//画第一条控制点之间的直线，用蓝色
	for (int i = 0;i+3<ControlPointTemp.size();i+=3){
		drawBezier(frame, ControlPointTemp[i], ControlPointTemp[i + 1], ControlPointTemp[i + 2], ControlPointTemp[i + 3], Scalar(0, 0, 255));

		//画连接控制点之间的直线，用蓝色
		line(frame, ControlPointTemp[i + 2], ControlPointTemp[i + 3], Scalar( 255,0, 0), 1);
		if (i+4<ControlPointTemp.size())
			line(frame, ControlPointTemp[i + 3], ControlPointTemp[i + 4], Scalar(255, 0, 0), 1);

	}
}





void QTandOpenCV::drawAllControlPointFlag(){
	vector<Point> ControlPointTemp(pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints);
	for (int i = 0; i < ControlPointTemp.size(); ++i){//缩放到显示大小,并对齐显示区域
		ControlPointTemp[i].x = ControlPointTemp[i].x *scale - hTranslation;
		ControlPointTemp[i].y = ControlPointTemp[i].y *scale - vTranslation;
	}
	ControlPointFlag = Mat(fmin(scaledHeight,height),fmin(scaledWidth, width), CV_32SC1,Scalar(-1));//int型，保存控制点的index，-1时不为控制点，每次更新
	for (int i = 0, ControlPointsSize = ControlPointTemp.size(); i < ControlPointsSize; ++i){//超过区域的地方会自动不画
		if (i == ControlPointsSize - 1 && ControlPointTemp[i] == ControlPointTemp[0])//最后一个控制点和第一个一样的话则标志位写0
			circle(ControlPointFlag, ControlPointTemp[i], FlagRadius, Scalar(0), -1);
		else
			circle(ControlPointFlag, ControlPointTemp[i], FlagRadius, Scalar(i), -1);
	}
}
void QTandOpenCV::drawControlPointFlag(int CurrentControlPoint,int x,int y){//直接获得的坐标
	circle(ControlPointFlag, Point(x, y), FlagRadius, Scalar(CurrentControlPoint), -1);
}
void QTandOpenCV::drawAllControlPoint(){
	cv::resize(frameOriginal, frame, Size(scaledWidth, scaledHeight), 0, 0, INTER_LINEAR);//和qt里面有冲突定义
	for (int i = 0; i < pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints.size(); ++i){
		Point temp = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[i];
		temp.x *= scale;
		temp.y *= scale;
		circle(frame,temp , ControlPointRadius, Scalar(0, 255, 0), -1);
	}
}
void QTandOpenCV::drawControlPoint(int x,int y){//输入的是缩放之后的坐标
	circle(frame, Point(x+hTranslation, y+vTranslation), ControlPointRadius, Scalar(0, 255, 0), -1);
}
void QTandOpenCV::updateAllControlPoint(int dx,int dy){//输入的是转换到原图尺寸的平移量
	int size = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints.size();
	for (int i = 0; i < size-1; ++i){
		if (i == 0 && size>1 && pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[0] == pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[size - 1]){//如果第一个控制点等于最后一个控制点
			pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[0].x += dx;
			pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[0].y += dy;
			pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[size - 1].x += dx;
			pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[size - 1].y += dy;
		}
		else{
			pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[i].x += dx;
			pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[i].y += dy;
		}
	}
	drawAllControlPointFlag();//重新画控制点的flag
}
void QTandOpenCV::updateControlPoint(int CurrentControlPoint, int dx,int dy){//输入的是缩放和滚动之后的坐标
	Point ControlPointTemp = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[CurrentControlPoint];
	ControlPointTemp.x = ControlPointTemp.x *scale - hTranslation;
	ControlPointTemp.y = ControlPointTemp.y *scale - vTranslation;
	circle(ControlPointFlag, ControlPointTemp, FlagRadius, Scalar(-1), -1);//抹掉原来的

	pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[CurrentControlPoint].x += dx;
	pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[CurrentControlPoint].y += dy;

	ControlPointTemp = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[CurrentControlPoint];
	ControlPointTemp.x = ControlPointTemp.x *scale - hTranslation;
	ControlPointTemp.y = ControlPointTemp.y *scale - vTranslation;
	circle(ControlPointFlag, ControlPointTemp, FlagRadius, Scalar(CurrentControlPoint), -1);//修改现在的
}
void QTandOpenCV::updateFirstControlPoint(int dx, int dy){
	Point ControlPointTemp = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[CurrentControlPoint];
	ControlPointTemp.x = ControlPointTemp.x *scale - hTranslation;
	ControlPointTemp.y = ControlPointTemp.y *scale - vTranslation;
	circle(ControlPointFlag, ControlPointTemp, FlagRadius, Scalar(-1), -1);//抹掉原来的

	int ControlPointsSize = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints.size();
	pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[0].x += dx;
	pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[0].y += dy;
	pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[ControlPointsSize-1].x += dx;
	pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[ControlPointsSize-1].y += dy;

	ControlPointTemp = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[CurrentControlPoint];
	ControlPointTemp.x = ControlPointTemp.x *scale - hTranslation;
	ControlPointTemp.y = ControlPointTemp.y *scale - vTranslation;
	circle(ControlPointFlag, ControlPointTemp, FlagRadius, Scalar(0), -1);//修改现在的
}
void QTandOpenCV::copyControlPoint(){//复制最邻近关键帧的控制点到当前关键帧
	int LastKeyFrame = CurrentFrame - 1;
	int NextKeyFrame = CurrentFrame + 1;
	while (pRoto->layers[CurrentLayer].KeyframeFlag[LastKeyFrame] == false){
		LastKeyFrame--;
	}
	while (pRoto->layers[CurrentLayer].KeyframeFlag[NextKeyFrame] == false){
		NextKeyFrame++;
	}
	int NearestKeyFrame;
	if (CurrentFrame-LastKeyFrame<=NextKeyFrame-CurrentFrame){
		NearestKeyFrame = LastKeyFrame;
	}
	else {
		NearestKeyFrame = NextKeyFrame;
	}
	pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints.assign(pRoto->layers[CurrentLayer].frames[NearestKeyFrame].ControlPoints.begin(),
		pRoto->layers[CurrentLayer].frames[NearestKeyFrame].ControlPoints.end());
	drawAllControlPointFlag();
}





void QTandOpenCV::setLinearFit(){
	if (ui.linearfit->isChecked()){
		InterpolationFlag = Linear;
	}
}
void QTandOpenCV::setParabolicFit(){
	if (ui.parabolicfit->isChecked()){
		InterpolationFlag = Parabolic;
	}
}
void QTandOpenCV::setCubicSplineFit(){
	if (ui.cubicsplinefit->isChecked()){
		InterpolationFlag = CubicSpline;
	}
}
void QTandOpenCV::setPolynomialFit(){
	if (ui.polynomialfit->isChecked()){
		InterpolationFlag = Polynomial;
		//如果视频还没开，则退出
		if (!video.isOpened() && files.size() == 0) return;
		int KeyframeSize = pRoto->layers[CurrentLayer].KeyframeIndex.size();
		if (KeyframeSize - 1>3)//n次多项式求解析方程需要n+1个点
			ui.npolynomialfit->setRange(3, KeyframeSize - 1);
	}
}
void QTandOpenCV::setPolynomialFitN(){
	if (InterpolationFlag == Polynomial){
		PolynomialFitN = ui.npolynomialfit->value();//n次多项式求解析方程需要n+1个点
	}
}



void QTandOpenCV::drawBezierInMask(Mat & frame, Point x0, Point x1, Point x2, Point x3, int type = 0){
	int BezierPointSize = bezierStepN(x0, x1, x2, x3);
	double step = 1.0 / (BezierPointSize + 1);

	double t = step;	//  只求中间的点

	frame.at<unsigned char>(x0.y, x0.x) = 255;
	frame.at<unsigned char>(x3.y, x3.x) = 255;
	int xPre = -1, yPre = -1;

	for (int i = 0; i<BezierPointSize; i++, t += step)
	{
		double coeff0 = (1 - t)*(1 - t)*(1 - t);
		double coeff1 = 3 * (1 - t)*(1 - t)*t;
		double coeff2 = 3 * (1 - t)*t*t;
		double coeff3 = t*t*t;
		double x = coeff0 * x0.x + coeff1 * x1.x + coeff2*x2.x + coeff3*x3.x;
		double y = coeff0 * x0.y + coeff1 * x1.y + coeff2*x2.y + coeff3*x3.y;
		int xNear = x + 0.5;//四舍五入，最近的点
		int yNear = y + 0.5;

		if (type == CV_AA){//反走样
			int xDown = x;//下界
			int yDown = y;
			int xUp = ceil(x);//上界
			int yUp = ceil(y);
			
			//距离最大为2^0.5
			float distUpL = pow(pow(yUp - y, 2) + pow(x - xDown, 2), 0.5);
			float distUpR = pow(pow(yUp - y, 2) + pow(xUp - x, 2), 0.5);
			float distDownL = pow(pow(y - yDown, 2) + pow(x - xDown, 2), 0.5);
			float distDownR = pow(pow(y - yDown, 2) + pow(xUp - x, 2), 0.5);

			float scaleUpL = 1 - 0.7*distUpL;//0.7 < 1/2^0.5
			float scaleUpR = 1 - 0.7*distUpR;
			float scaleDownL = 1 - 0.7*distDownL;
			float scaleDownR = 1 - 0.7*distDownR;

			float tempUpL = frame.at<unsigned char>(yUp, xDown);
			float tempUpR = frame.at<unsigned char>(yUp, xUp);
			float tempDownL = frame.at<unsigned char>(yDown, xDown);
			float tempDownR = frame.at<unsigned char>(yDown, xUp);

			

			//if (xPre == xNear && yPre == yNear){
				if (xyTest(xDown, yUp, frame.cols, frame.rows)) frame.at<unsigned char>(yUp, xDown) = max(tempUpL, scaleUpL * 255);
				if (xyTest(xUp, yUp, frame.cols, frame.rows)) frame.at<unsigned char>(yUp, xUp) = max(tempUpR, scaleUpR * 255);
				if (xyTest(xUp, yDown, frame.cols, frame.rows)) frame.at<unsigned char>(yDown, xDown) = max(tempDownL, scaleDownL * 255);
				if (xyTest(xUp, yUp, frame.cols, frame.rows)) frame.at<unsigned char>(yDown, xUp) = max(tempDownR, scaleDownR * 255);
			//}
			/*else{
				if (xyTest(xDown, yUp, frame.cols, frame.rows)) frame.at<unsigned char>(yUp, xDown) = scaleUpL * 255;
				if (xyTest(xUp, yUp, frame.cols, frame.rows)) frame.at<unsigned char>(yUp, xUp) = scaleUpR * 255;
				if (xyTest(xUp, yDown, frame.cols, frame.rows)) frame.at<unsigned char>(yDown, xDown) = scaleDownL * 255;
				if (xyTest(xUp, yUp, frame.cols, frame.rows)) frame.at<unsigned char>(yDown, xUp) = scaleDownR * 255;
				xPre = xNear;
				yPre = yNear;
			}*/

			
		}
		else{//取最近点
			double distance = pow(pow(x - xNear, 2) + pow(y - yNear, 2), 0.5);
			if (xyTest(xNear, yNear, frame.cols, frame.rows)) frame.at<unsigned char>(yNear, xNear) = 255;// *(1 - distance);
		}
	}
}

Mat QTandOpenCV::getMask(int iCurrentLayer, int iCurrentFrame){
	Mat mask(Height, Width, CV_8UC1, Scalar(0));
	vector<Point> ControlPointTemp(pRoto->layers[iCurrentLayer].frames[iCurrentFrame].ControlPoints);
	for (int i = 0; i + 3<ControlPointTemp.size(); i += 3){
		drawBezierInMask(mask,ControlPointTemp[i], ControlPointTemp[i + 1], ControlPointTemp[i + 2], ControlPointTemp[i + 3]);
	}

	MemStorage Storage;
	vector<vector<Point> > contours;
	// 查找所有轮廓 
	//void findContours(InputOutputArray image, OutputArrayOfArrays contours,int mode, int method, Point offset = Point());
	//findContours(mask, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);//保存轮廓直线
	findContours(mask, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);//保存轮廓点


	// 填充所有轮廓,lineType,8：八连通线型，4：bresenham，4-connected
	//void drawContours(InputOutputArray image, InputArrayOfArrays contours, int contourIdx, const Scalar& color, int thickness = 1, int lineType = 8,	InputArray hierarchy=noArray(),	int maxLevel=INT_MAX, Point offset=Point() );
	//drawContours(mask, contours, -1, Scalar(255), CV_FILLED,CV_AA);//采用反走样方式来画直线
	drawContours(mask, contours, -1, Scalar(255), CV_FILLED);//采用反走样方式来画直线
	

	//腐蚀erodes the image (applies the local minimum operator)
	//void erode( InputArray src, OutputArray dst, InputArray kernel, Point anchor=Point(-1,-1), int iterations=1, int borderType=BORDER_CONSTANT, const Scalar& borderValue=morphologyDefaultBorderValue() );
	//erode(mask, mask, NULL);

	//膨胀dilates the image (applies the local maximum operator)
	//void dilate( InputArray src, OutputArray dst, InputArray kernel, Point anchor=Point(-1,-1), int iterations=1, int borderType=BORDER_CONSTANT, const Scalar& borderValue=morphologyDefaultBorderValue() );
	//dilate(mask,mask,NULL);


	//重画反走样轮廓
	for (int i = 0; i + 3<ControlPointTemp.size(); i += 3){
		drawBezierInMask(mask, ControlPointTemp[i], ControlPointTemp[i + 1], ControlPointTemp[i + 2], ControlPointTemp[i + 3],CV_AA);
	}
	//blur(mask,mask,Size(2,2),Point(-1,-1));
	//GaussianBlur(mask,mask,Size(3,3),0,0,BORDER_DEFAULT);

	return mask;
}

Mat QTandOpenCV::getMaskBGRA(int iCurrentLayer, int iCurrentFrame){
	Mat mask = getMask(iCurrentLayer,  iCurrentFrame);
	Mat blend(Height, Width, CV_8UC4, Scalar(0));
	Mat temp;
	if (inputType == Video){
		video.set(CV_CAP_PROP_POS_FRAMES, iCurrentFrame);
		video >> temp;
	}
	else if (inputType == Sequence){
		temp = imread(files[iCurrentFrame]);
	}
	vector<Mat> channels;
	split(temp,channels);
	channels.push_back(mask);
	merge(channels,blend);
	return blend;
}

void QTandOpenCV::getResult(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;

	//初始化非关键帧的控制点
	for (int i = 0; i < pRoto->layers.size(); ++i){
		for (int j = 0; j < pRoto->layers[i].frames.size(); ++j){
			if (!pRoto->layers[i].KeyframeFlag[j]){
				pRoto->layers[i].frames[j].ControlPoints.assign(pRoto->layers[i].frames[0].ControlPoints.begin(),
					pRoto->layers[i].frames[0].ControlPoints.end());
			}
		}
		pRoto->layers[i].GetResultDone = true;
	}
	if (InterpolationFlag == Linear){
		for (int i = 0; i < pRoto->layers.size(); ++i){//layer
			int ControlPointNum = pRoto->layers[i].frames[0].ControlPoints.size();
			for (int n = 0; n < ControlPointNum; ++n){//control point
				for (int j = 0; j < pRoto->layers[i].KeyframeIndex.size() - 1; ++j){//key frame 
					int keyframe1 = pRoto->layers[i].KeyframeIndex[j];
					int keyframe2 = pRoto->layers[i].KeyframeIndex[j+1];
					int x1 = pRoto->layers[i].frames[keyframe1].ControlPoints[n].x;
					int y1 = pRoto->layers[i].frames[keyframe1].ControlPoints[n].y;
					int x2 = pRoto->layers[i].frames[keyframe2].ControlPoints[n].x;
					int y2 = pRoto->layers[i].frames[keyframe2].ControlPoints[n].y;
					float gradient_x = (float(x2 - x1)) / (keyframe2 - keyframe1);
					float b_x = x1 - keyframe1*gradient_x;
					float gradient_y = (float(y2 - y1)) / (keyframe2 - keyframe1);
					float b_y = y1 - keyframe1*gradient_y;
					for (int k = keyframe1 + 1; k < keyframe2; ++k){//frame between key frame
						pRoto->layers[i].frames[k].ControlPoints[n].x = k*gradient_x + b_x;
						pRoto->layers[i].frames[k].ControlPoints[n].y = k*gradient_y + b_y;
					}
				}
			}
			pRoto->layers[i].GetResultDone = true;
		}
	}
	else if (InterpolationFlag == Parabolic){//不平滑，所以不采用
		for (int i = 0; i < pRoto->layers.size(); ++i){//layer
			int ControlPointNum = pRoto->layers[i].frames[0].ControlPoints.size();
			int KeyframeSize = pRoto->layers[i].KeyframeIndex.size();
			for (int n = 0; n < ControlPointNum; ++n){//control point
				for (int j = 0; j < KeyframeSize - 2; j+=2){//key frame 
					int t1 = pRoto->layers[i].KeyframeIndex[j];
					int t2 = pRoto->layers[i].KeyframeIndex[j + 1];
					int t3 = pRoto->layers[i].KeyframeIndex[j + 2];
					double t12 = t1*t1;
					double t22 = t2*t2;
					double t32 = t3*t3;
					double fenmu = (t1 - t3)*(t1 - t2)*(t2 - t3);


					double x1 = pRoto->layers[i].frames[t1].ControlPoints[n].x;
					double x2 = pRoto->layers[i].frames[t2].ControlPoints[n].x;
					double x3 = pRoto->layers[i].frames[t3].ControlPoints[n].x;
					double a_x = (t1 * x3 + t2 * x1 + t3 * x2 - t1 * x2 - t2 * x3 - t3 * x1) / fenmu;
					double b_x = (t12 * x2 + t22 * x3 + t32 * x1 - t12 * x3 - t22 * x1 - t32 * x2) / fenmu;
					double c_x = x1 - a_x*t1*t1 - b_x*t1;


					double y1 = pRoto->layers[i].frames[t1].ControlPoints[n].y;
					double y2 = pRoto->layers[i].frames[t2].ControlPoints[n].y;
					double y3 = pRoto->layers[i].frames[t3].ControlPoints[n].y;
					double a_y = (t1 * y3 + t2 * y1 + t3 * y2 - t1 * y2 - t2 * y3 - t3 * y1) / fenmu;
					double b_y = (t12 * y2 + t22 * y3 + t32 * y1 - t12 * y3 - t22 * y1 - t32 * y2) / fenmu;
					double c_y = y1 - a_y*t1*t1 - b_y*t1;


					for (int k = t1 + 1; k < t2; ++k){//frame between key frame t1-t2
						int x = a_x*k*k + b_x*k + c_x;
						int y = a_y*k*k + b_y*k + c_y;
						pRoto->layers[i].frames[k].ControlPoints[n].x = x;
						pRoto->layers[i].frames[k].ControlPoints[n].y = y;
					}
					for (int k = t2 + 1; k < t3; ++k){//frame between key frame t2-t3
						int x = a_x*k*k + b_x*k + c_x;
						int y = a_y*k*k + b_y*k + c_y;
						pRoto->layers[i].frames[k].ControlPoints[n].x = x;
						pRoto->layers[i].frames[k].ControlPoints[n].y = y;
					}
				}
			}
			pRoto->layers[i].GetResultDone = true;
		}
	}
	else if (InterpolationFlag == CubicSpline){//自由边界
		for (int i = 0; i < pRoto->layers.size(); ++i){//layer
			int ControlPointNum = pRoto->layers[i].frames[0].ControlPoints.size();
			int KeyframeSize = pRoto->layers[i].KeyframeIndex.size();
			for (int n = 0; n < ControlPointNum; ++n){//control point
				vector<float> h(KeyframeSize - 1, 0);
				vector<float> a(KeyframeSize, 0);
				vector<float> b(KeyframeSize, 1);
				vector<float> c(KeyframeSize - 1, 0);
				vector<float> d_x(KeyframeSize, 0);
				vector<float> d_y(KeyframeSize, 0);
				vector<float> m_x(KeyframeSize, 0);
				vector<float> m_y(KeyframeSize, 0);
				for (int j = 0; j < KeyframeSize - 1; ++j){
					h[j] = pRoto->layers[i].KeyframeIndex[j + 1] - pRoto->layers[i].KeyframeIndex[j];
				}
				for (int j = 1; j < KeyframeSize - 1;++j){
					b[j] = 2 * (pRoto->layers[i].KeyframeIndex[j + 1] - pRoto->layers[i].KeyframeIndex[j - 1]);
				}
				for (int j = 1; j < KeyframeSize - 1; ++j){
					a[j] = pRoto->layers[i].KeyframeIndex[j] - pRoto->layers[i].KeyframeIndex[j - 1];
				}
				for (int j = 1; j < KeyframeSize - 1; ++j){
					c[j] = pRoto->layers[i].KeyframeIndex[j + 1] - pRoto->layers[i].KeyframeIndex[j];
				}
				for (int j = 1; j < KeyframeSize - 1; ++j){
					d_x[j] = 6 * ((pRoto->layers[i].frames[pRoto->layers[i].KeyframeIndex[j + 1]].ControlPoints[n].x - pRoto->layers[i].frames[pRoto->layers[i].KeyframeIndex[j]].ControlPoints[n].x) / h[j] -
						(pRoto->layers[i].frames[pRoto->layers[i].KeyframeIndex[j]].ControlPoints[n].x - pRoto->layers[i].frames[pRoto->layers[i].KeyframeIndex[j-1]].ControlPoints[n].x) / h[j - 1]);
				}
				for (int j = 1; j < KeyframeSize - 1; ++j){
					d_y[j] = 6 * ((pRoto->layers[i].frames[pRoto->layers[i].KeyframeIndex[j + 1]].ControlPoints[n].y - pRoto->layers[i].frames[pRoto->layers[i].KeyframeIndex[j]].ControlPoints[n].y) / h[j] -
						(pRoto->layers[i].frames[pRoto->layers[i].KeyframeIndex[j]].ControlPoints[n].y - pRoto->layers[i].frames[pRoto->layers[i].KeyframeIndex[j - 1]].ControlPoints[n].y) / h[j - 1]);
				}
				tdma(m_x, a, b, c, d_x);
				tdma(m_y, a, b, c, d_y);


				for (int j = 0; j < pRoto->layers[i].KeyframeIndex.size() - 1; ++j){//key frame 
					int keyframe1 = pRoto->layers[i].KeyframeIndex[j];
					int keyframe2 = pRoto->layers[i].KeyframeIndex[j + 1];

					float ai_x = pRoto->layers[i].frames[keyframe1].ControlPoints[n].x;
					float bi_x = (pRoto->layers[i].frames[keyframe2].ControlPoints[n].x - pRoto->layers[i].frames[keyframe1].ControlPoints[n].x) / h[j] - h[j] * m_x[j] / 2 - h[j] * (m_x[j + 1] - m_x[j]) / 6;
					float ci_x = m_x[j] / 2;
					float di_x = (m_x[j + 1] - m_x[j]) / (6 * h[j]);

					float ai_y = pRoto->layers[i].frames[keyframe1].ControlPoints[n].y;
					float bi_y = (pRoto->layers[i].frames[keyframe2].ControlPoints[n].y - pRoto->layers[i].frames[keyframe1].ControlPoints[n].y) / h[j] - h[j] * m_y[j] / 2 - h[j] * (m_y[j + 1] - m_y[j]) / 6;
					float ci_y = m_y[j] / 2;
					float di_y = (m_y[j + 1] - m_y[j]) / (6 * h[j]);

					for (int k = keyframe1 + 1; k < keyframe2; ++k){//frame between key frame
						float x_xi = k - keyframe1;
						float x_xi_2 = x_xi * x_xi;
						float x_xi_3 = x_xi_2 * x_xi;
						pRoto->layers[i].frames[k].ControlPoints[n].x = ai_x + bi_x*x_xi + ci_x*x_xi_2 + di_x*x_xi_3;
						pRoto->layers[i].frames[k].ControlPoints[n].y = ai_y + bi_y*x_xi + ci_y*x_xi_2 + di_y*x_xi_3;
					}
				}
			}
			pRoto->layers[i].GetResultDone = true;
		}
	}
	else if (InterpolationFlag == Polynomial){//效果不好
		
		for (int i = 0; i < pRoto->layers.size(); ++i){//layer
			int ControlPointNum = pRoto->layers[i].frames[0].ControlPoints.size();
			int KeyframeSize = pRoto->layers[i].KeyframeIndex.size();
			vector<float> t(KeyframeSize);
			vector<float> x(KeyframeSize);
			vector<float> y(KeyframeSize);
			vector<float> ax(PolynomialFitN+1,0);//初始化
			vector<float> ay(PolynomialFitN+1,0);//初始化
			vector<float> dtx(3);
			vector<float> dty(3);

			for (int n = 0; n < ControlPointNum; ++n){//control point
				for (int j = 0; j < KeyframeSize; ++j){
					int KeyframeIndex = pRoto->layers[i].KeyframeIndex[j];
					t[j] = KeyframeIndex;
					x[j] = pRoto->layers[i].frames[KeyframeIndex].ControlPoints[n].x;
					y[j] = pRoto->layers[i].frames[KeyframeIndex].ControlPoints[n].y;
				}

				curvePolyFit(t, x, KeyframeSize, ax, PolynomialFitN + 1, dtx);
				curvePolyFit(t, y, KeyframeSize, ay, PolynomialFitN + 1, dty);


				for (int j = 0; j < pRoto->layers[i].KeyframeIndex.size() - 1; ++j){//key frame 
					int keyframe1 = pRoto->layers[i].KeyframeIndex[j];
					int keyframe2 = pRoto->layers[i].KeyframeIndex[j + 1];
					for (int k = keyframe1 + 1; k < keyframe2; ++k){//frame between key frame
						pRoto->layers[i].frames[k].ControlPoints[n].x = curvePolyFitValue(k, ax, PolynomialFitN + 1, KeyframeSize);
						pRoto->layers[i].frames[k].ControlPoints[n].y = curvePolyFitValue(k, ay, PolynomialFitN + 1, KeyframeSize);
					}
				}
			}
			pRoto->layers[i].GetResultDone = true;
		}

		//void polyfit(const Mat& srcx, const Mat& srcy, Mat& dst, int order);
		//polyfit();

	}

	
	drawAllControlPointFlag();
	drawAllControlPoint();
	drawCurve();
	showFrame(frame);
	QMessageBox::warning(this, tr("get result"), tr("get result done!"));
}

void QTandOpenCV::updateResultLinear(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	if (!pRoto->layers[CurrentLayer].KeyframeFlag[CurrentFrame]) return;
	int temp = 0;
	while (pRoto->layers[CurrentLayer].KeyframeIndex[temp] != CurrentFrame){
		temp++;
	}
	pRoto->layers[CurrentLayer].CurrentKeyframeIndex = temp;
	//用线性插值来暂时更新当前被修改的关键帧附近的帧
	int ControlPointNum = pRoto->layers[CurrentLayer].frames[0].ControlPoints.size();
	for (int n = 0; n < ControlPointNum; ++n){//control point
		int KeyFrameSize = pRoto->layers[CurrentLayer].KeyframeIndex.size();
		for (int j = max(0, pRoto->layers[CurrentLayer].CurrentKeyframeIndex - 1); j < max(KeyFrameSize, pRoto->layers[CurrentLayer].CurrentKeyframeIndex + 1) - 1; ++j){//key frame 
			int keyframe1 = pRoto->layers[CurrentLayer].KeyframeIndex[j];
			int keyframe2 = pRoto->layers[CurrentLayer].KeyframeIndex[j + 1];
			int x1 = pRoto->layers[CurrentLayer].frames[keyframe1].ControlPoints[n].x;
			int y1 = pRoto->layers[CurrentLayer].frames[keyframe1].ControlPoints[n].y;
			int x2 = pRoto->layers[CurrentLayer].frames[keyframe2].ControlPoints[n].x;
			int y2 = pRoto->layers[CurrentLayer].frames[keyframe2].ControlPoints[n].y;
			float gradient_x = (float(x2 - x1)) / (keyframe2 - keyframe1);
			float b_x = x1 - keyframe1*gradient_x;
			float gradient_y = (float(y2 - y1)) / (keyframe2 - keyframe1);
			float b_y = y1 - keyframe1*gradient_y;
			for (int k = keyframe1 + 1; k < keyframe2; ++k){//frame between key frame
				pRoto->layers[CurrentLayer].frames[k].ControlPoints[n].x = k*gradient_x + b_x;
				pRoto->layers[CurrentLayer].frames[k].ControlPoints[n].y = k*gradient_y + b_y;
			}
		}
	}
	//pRoto->layers[CurrentLayer].GetResultDone = true;
}

void QTandOpenCV::updateResult(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	if (!pRoto->layers[CurrentLayer].KeyframeFlag[CurrentFrame]) return;
	int temp = 0;
	while (pRoto->layers[CurrentLayer].KeyframeIndex[temp] != CurrentFrame){
		temp++;
	}
	pRoto->layers[CurrentLayer].CurrentKeyframeIndex = temp;
	int ControlPointNum = pRoto->layers[CurrentLayer].frames[0].ControlPoints.size();
	int KeyframeSize = pRoto->layers[CurrentLayer].KeyframeIndex.size();
	if (InterpolationFlag == Linear){
		int ControlPointNum = pRoto->layers[CurrentLayer].frames[0].ControlPoints.size();
		for (int n = 0; n < ControlPointNum; ++n){//control point
			for (int j = max(0, pRoto->layers[CurrentLayer].CurrentKeyframeIndex - 1); j < max(KeyframeSize, pRoto->layers[CurrentLayer].CurrentKeyframeIndex + 1) - 1; ++j){//key frame 
				int keyframe1 = pRoto->layers[CurrentLayer].KeyframeIndex[j];
				int keyframe2 = pRoto->layers[CurrentLayer].KeyframeIndex[j + 1];
				int x1 = pRoto->layers[CurrentLayer].frames[keyframe1].ControlPoints[n].x;
				int y1 = pRoto->layers[CurrentLayer].frames[keyframe1].ControlPoints[n].y;
				int x2 = pRoto->layers[CurrentLayer].frames[keyframe2].ControlPoints[n].x;
				int y2 = pRoto->layers[CurrentLayer].frames[keyframe2].ControlPoints[n].y;
				float gradient_x = (float(x2 - x1)) / (keyframe2 - keyframe1);
				float b_x = x1 - keyframe1*gradient_x;
				float gradient_y = (float(y2 - y1)) / (keyframe2 - keyframe1);
				float b_y = y1 - keyframe1*gradient_y;
				for (int k = keyframe1 + 1; k < keyframe2; ++k){//frame between key frame
					pRoto->layers[CurrentLayer].frames[k].ControlPoints[n].x = k*gradient_x + b_x;
					pRoto->layers[CurrentLayer].frames[k].ControlPoints[n].y = k*gradient_y + b_y;
				}
			}
		}
	}
	else if (InterpolationFlag == CubicSpline){//自由边界
		for (int n = 0; n < ControlPointNum; ++n){//control point
			vector<float> h(KeyframeSize - 1, 0);
			vector<float> a(KeyframeSize, 0);
			vector<float> b(KeyframeSize, 1);
			vector<float> c(KeyframeSize - 1, 0);
			vector<float> d_x(KeyframeSize, 0);
			vector<float> d_y(KeyframeSize, 0);
			vector<float> m_x(KeyframeSize, 0);
			vector<float> m_y(KeyframeSize, 0);
			for (int j = 0; j < KeyframeSize - 1; ++j){
				h[j] = pRoto->layers[CurrentLayer].KeyframeIndex[j + 1] - pRoto->layers[CurrentLayer].KeyframeIndex[j];
			}
			for (int j = 1; j < KeyframeSize - 1; ++j){
				b[j] = 2 * (pRoto->layers[CurrentLayer].KeyframeIndex[j + 1] - pRoto->layers[CurrentLayer].KeyframeIndex[j - 1]);
			}
			for (int j = 1; j < KeyframeSize - 1; ++j){
				a[j] = pRoto->layers[CurrentLayer].KeyframeIndex[j] - pRoto->layers[CurrentLayer].KeyframeIndex[j - 1];
			}
			for (int j = 1; j < KeyframeSize - 1; ++j){
				c[j] = pRoto->layers[CurrentLayer].KeyframeIndex[j + 1] - pRoto->layers[CurrentLayer].KeyframeIndex[j];
			}
			for (int j = 1; j < KeyframeSize - 1; ++j){
				d_x[j] = 6 * ((pRoto->layers[CurrentLayer].frames[pRoto->layers[CurrentLayer].KeyframeIndex[j + 1]].ControlPoints[n].x - pRoto->layers[CurrentLayer].frames[pRoto->layers[CurrentLayer].KeyframeIndex[j]].ControlPoints[n].x) / h[j] -
					(pRoto->layers[CurrentLayer].frames[pRoto->layers[CurrentLayer].KeyframeIndex[j]].ControlPoints[n].x - pRoto->layers[CurrentLayer].frames[pRoto->layers[CurrentLayer].KeyframeIndex[j - 1]].ControlPoints[n].x) / h[j - 1]);
			}
			for (int j = 1; j < KeyframeSize - 1; ++j){
				d_y[j] = 6 * ((pRoto->layers[CurrentLayer].frames[pRoto->layers[CurrentLayer].KeyframeIndex[j + 1]].ControlPoints[n].y - pRoto->layers[CurrentLayer].frames[pRoto->layers[CurrentLayer].KeyframeIndex[j]].ControlPoints[n].y) / h[j] -
					(pRoto->layers[CurrentLayer].frames[pRoto->layers[CurrentLayer].KeyframeIndex[j]].ControlPoints[n].y - pRoto->layers[CurrentLayer].frames[pRoto->layers[CurrentLayer].KeyframeIndex[j - 1]].ControlPoints[n].y) / h[j - 1]);
			}
			tdma(m_x, a, b, c, d_x);
			tdma(m_y, a, b, c, d_y);


			for (int j = max(0, pRoto->layers[CurrentLayer].CurrentKeyframeIndex - 1); j < max(KeyframeSize, pRoto->layers[CurrentLayer].CurrentKeyframeIndex + 1) - 1; ++j){//key frame 
				int keyframe1 = pRoto->layers[CurrentLayer].KeyframeIndex[j];
				int keyframe2 = pRoto->layers[CurrentLayer].KeyframeIndex[j + 1];

				float ai_x = pRoto->layers[CurrentLayer].frames[keyframe1].ControlPoints[n].x;
				float bi_x = (pRoto->layers[CurrentLayer].frames[keyframe2].ControlPoints[n].x - pRoto->layers[CurrentLayer].frames[keyframe1].ControlPoints[n].x) / h[j] - h[j] * m_x[j] / 2 - h[j] * (m_x[j + 1] - m_x[j]) / 6;
				float ci_x = m_x[j] / 2;
				float di_x = (m_x[j + 1] - m_x[j]) / (6 * h[j]);

				float ai_y = pRoto->layers[CurrentLayer].frames[keyframe1].ControlPoints[n].y;
				float bi_y = (pRoto->layers[CurrentLayer].frames[keyframe2].ControlPoints[n].y - pRoto->layers[CurrentLayer].frames[keyframe1].ControlPoints[n].y) / h[j] - h[j] * m_y[j] / 2 - h[j] * (m_y[j + 1] - m_y[j]) / 6;
				float ci_y = m_y[j] / 2;
				float di_y = (m_y[j + 1] - m_y[j]) / (6 * h[j]);

				for (int k = keyframe1 + 1; k < keyframe2; ++k){//frame between key frame
					float x_xi = k - keyframe1;
					float x_xi_2 = x_xi * x_xi;
					float x_xi_3 = x_xi_2 * x_xi;
					pRoto->layers[CurrentLayer].frames[k].ControlPoints[n].x = ai_x + bi_x*x_xi + ci_x*x_xi_2 + di_x*x_xi_3;
					pRoto->layers[CurrentLayer].frames[k].ControlPoints[n].y = ai_y + bi_y*x_xi + ci_y*x_xi_2 + di_y*x_xi_3;
				}
			}
		}
	}
	else if (InterpolationFlag == Parabolic){//不平滑，所以不采用
	}
	else if (InterpolationFlag == Polynomial){//效果不好
	}
	//pRoto->layers[CurrentLayer].GetResultDone = true;
}


void QTandOpenCV::createFile(string Path){
	QDir *temp = new QDir;
	if (!temp->exists(QString::fromStdString(Path)))
		temp->mkdir(QString::fromStdString(Path));

	/*if ( temp->exists( QString::fromStdString(Path) ) )
		QMessageBox::warning(this, tr("create"), tr("exist!"));
	else{
		bool ok = temp->mkdir( QString::fromStdString(Path) );
		if (ok)
			QMessageBox::warning(this, tr("create"), tr("success!"));
	}*/
}
void QTandOpenCV::deleteFile(string Path){
	QDir *temp = new QDir;
	if (temp->exists(QString::fromStdString(Path)))
		temp->rmdir(QString::fromStdString(Path));
}
void QTandOpenCV::saveResultText(string Path){
	ofstream os;

	//std::locale::global(std::locale(""));//将全局区域设为操作系统默认区域//可以顺利打开文件了
	//os.open(Path.c_str(), ofstream::out);
	//std::locale::global(std::locale("C"));//还原全局区域设定
	
	//setlocale(LC_ALL, "Chinese-simplified");//设置中文环境
	//os.open(Path, ofstream::out);
	//setlocale(LC_ALL, "C");//还原


	os.open((const char *)QString::fromStdString(Path).toLocal8Bit());
	if (!os) {
		ui.label->setText(QString::fromStdString("Can not open file:\n"+Path));
		return;
	}
	os << "<FrameNum>\n" << pRoto->FrameNum << "\n</FrameNum>\n" << endl;
	os << "<Height>\n" << pRoto->Height << "\n</Height>\n" << endl;
	os << "<Width>\n" << pRoto->Width << "\n</Width>\n" << endl;
	os << "<Layers>\n" << pRoto->layers.size() << "\n</Layers>\n" << endl;
	for (int i = 0; i < pRoto->layers.size();++i){
		os << "<Layer>\n" << i << "\n</Layer>\n" << endl;
		os << "<KeyframeNum>\n" << pRoto->layers[i].KeyframeIndex.size() << "\n</KeyframeNum>\n" << endl;
		for (int j = 0; j < pRoto->layers[i].KeyframeIndex.size();++j){
			os << " " << pRoto->layers[i].KeyframeIndex[j];
		}
		os << "\n" << endl;
		os << "<ControlPointNum>\n" << pRoto->layers[i].frames[0].ControlPoints.size() << "\n</ControlPointNum>\n" << endl;
		for (int j = 0; j < pRoto->FrameNum;++j){
			for (int k = 0; k < pRoto->layers[i].frames[0].ControlPoints.size();++k){
				os << " " << pRoto->layers[i].frames[j].ControlPoints[k].x << " " << pRoto->layers[i].frames[j].ControlPoints[k].y;
			}
			os << endl;
		}
		os << endl;
	}
	os.close();
}
void QTandOpenCV::readResultText(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;



	QString resultPath = QFileDialog::getOpenFileName(this, QString::fromLocal8Bit("请选择相应视频的已有结果文件："), "D:/", tr("text files(*.txt)"));
	ifstream is;
	is.open((const char *)resultPath.toLocal8Bit());
	if (!is) {
		ui.label->setText(QString::fromStdString("Can not open file:\n" + resultPath.toStdString()));
		return;
	}
	string tempString;
	int frameNum, height, width;
	is >> tempString >> frameNum >> tempString;
	is >> tempString >> height >> tempString;
	is >> tempString >> width >> tempString;
	if (pRoto->FrameNum != frameNum || pRoto->Height != height || pRoto->Width != width){
		QMessageBox::warning(this, tr("load result"), tr("not match!"));
		return;
	}
	int layerSize = 1;
	is >> tempString >> layerSize >> tempString;
	pRoto->layers.clear();//清理一下先
	for (int i = 0; i < layerSize;++i){
		addLayer();
		int layerIndex;
		is >> tempString >> layerIndex >> tempString;
		pRoto->layers[layerIndex].FirstKeyframeDone = true;
		pRoto->layers[layerIndex].GetResultDone = true;
		int keyframeNum;
		is >> tempString >> keyframeNum >> tempString;
		pRoto->layers[layerIndex].KeyframeIndex.clear();//清理一下先
		for (int j = 0; j < keyframeNum;++j){
			int keyframeIndex;
			is >>keyframeIndex ;
			pRoto->layers[layerIndex].KeyframeIndex.push_back(keyframeIndex);
			pRoto->layers[layerIndex].KeyframeFlag[keyframeIndex] = true;
		}
		int controlPointNum;
		is >> tempString >> controlPointNum >> tempString;
		for (int j = 0; j < pRoto->FrameNum; ++j){
			for (int k = 0; k < controlPointNum;++k){
				int x, y;
				is >> x >> y;
				pRoto->layers[i].frames[j].ControlPoints.push_back(Point(x, y));
			}
		}
	}
	QMessageBox::warning(this, tr("load result"), tr("load result done!"));
	drawAllControlPointFlag();
	drawAllControlPoint();
	drawCurve();
	showFrame(frame);
	showKeyFrameMark();

}
void QTandOpenCV::saveResult(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;






	string resultPath;//默认的保存路径
	if (inputType == Video){
		if (videoPath[videoPath.size() - 1] == '/'){
			resultPath = videoPath + videoName;
		}
		else{
			resultPath = videoPath +"/"+ videoName;
		}
	}
	else if(inputType==Sequence) {
		resultPath = videoPath + '_' + videoName;
	}
	createFile(resultPath);
	


	//选择已有的路径
	QString file_full = QFileDialog::getExistingDirectory(this, QString::fromLocal8Bit("请选择保存路径："), QString::fromStdString(resultPath));
	if (file_full.isEmpty()) return;
	if (QString::fromStdString(resultPath) != file_full) {
		deleteFile(resultPath);
		resultPath = file_full.toStdString();
	}

	

	//设置进度条
	int TotalFrame = pRoto->layers.size()*pRoto->FrameNum;
	int WorkCount = 0;
	QProgressDialog ProgressDialog(tr("save result"), tr("cancel"), 0, TotalFrame, this);
	ProgressDialog.setWindowTitle(tr("please wait..."));
	ProgressDialog.setWindowModality(Qt::WindowModal);
	ProgressDialog.show();





	for (int i = 0; i < pRoto->layers.size();++i){
		if (pRoto->layers[i].GetResultDone != true) {
			QMessageBox::warning(this, tr("save result"), tr("fail! please get result first!"));
			return;
		}

		string layerPath = resultPath + '/' + QString::number(i).toStdString() ;
		createFile(layerPath);
		for (int j = 0; j < pRoto->layers[i].frames.size();++j){
			char number[10] = {0};
			sprintf(number, "%08d", j);
			string name = number;
			string maskPath = layerPath + '/' + name + ".png";
			Mat mask = getMaskBGRA(i,j);
			QString temp = QString::fromStdString(maskPath);
			IplImage ipl_mask(mask);
			int params[] = { CV_IMWRITE_PNG_COMPRESSION, 9 };
			cvSaveImage((const char *)temp.toLocal8Bit(), &ipl_mask,params);//保存四通道图片


			//更新进度条
			ProgressDialog.setValue(++WorkCount);
			QCoreApplication::processEvents();
			if (ProgressDialog.wasCanceled()){
				return;
			}
				


			//string maskPath = layerPath + '/' + QString::number(j).toStdString() + ".png";
			//Mat mask = getMask(i,j);
			//imwrite(maskPath,mask);//不能保存中文路径
		}
	}

	//关闭进度条
	ProgressDialog.close();


	string resultPathText = resultPath + "/TextResult.txt";
	cout << resultPathText << endl;
	saveResultText(resultPathText);

	QMessageBox::warning(this, tr("save result"), tr("save result done!"));
}





void QTandOpenCV::nextKeyframe(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;

	int CurrentKeyFrame = pRoto->layers[CurrentLayer].KeyframeIndex[pRoto->layers[CurrentLayer].CurrentKeyframeIndex];
	if (CurrentKeyFrame <= CurrentFrame){
		while (pRoto->layers[CurrentLayer].KeyframeIndex[pRoto->layers[CurrentLayer].CurrentKeyframeIndex] <= CurrentFrame
			&& pRoto->layers[CurrentLayer].CurrentKeyframeIndex < pRoto->layers[CurrentLayer].KeyframeIndex.size() - 1){
			pRoto->layers[CurrentLayer].CurrentKeyframeIndex++;
		}
	}
	else{
		while (pRoto->layers[CurrentLayer].KeyframeIndex[pRoto->layers[CurrentLayer].CurrentKeyframeIndex] > CurrentFrame
			&& pRoto->layers[CurrentLayer].CurrentKeyframeIndex >0){
			pRoto->layers[CurrentLayer].CurrentKeyframeIndex--;
		}
		pRoto->layers[CurrentLayer].CurrentKeyframeIndex++;
	}
	CurrentFrame = pRoto->layers[CurrentLayer].KeyframeIndex[pRoto->layers[CurrentLayer].CurrentKeyframeIndex];
	ui.horizontalSlider->setSliderPosition(CurrentFrame);
}
void QTandOpenCV::lastKeyframe(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	int CurrentKeyFrame = pRoto->layers[CurrentLayer].KeyframeIndex[pRoto->layers[CurrentLayer].CurrentKeyframeIndex];
	if (CurrentKeyFrame >= CurrentFrame){
		while (pRoto->layers[CurrentLayer].KeyframeIndex[pRoto->layers[CurrentLayer].CurrentKeyframeIndex] >= CurrentFrame
			&& pRoto->layers[CurrentLayer].CurrentKeyframeIndex >0){
			pRoto->layers[CurrentLayer].CurrentKeyframeIndex--;
		}
	}
	else{
		while (pRoto->layers[CurrentLayer].KeyframeIndex[pRoto->layers[CurrentLayer].CurrentKeyframeIndex] < CurrentFrame
			&& pRoto->layers[CurrentLayer].CurrentKeyframeIndex < pRoto->layers[CurrentLayer].KeyframeIndex.size() - 1){
			pRoto->layers[CurrentLayer].CurrentKeyframeIndex++;
		}
		pRoto->layers[CurrentLayer].CurrentKeyframeIndex--;
	}
	CurrentFrame = pRoto->layers[CurrentLayer].KeyframeIndex[pRoto->layers[CurrentLayer].CurrentKeyframeIndex];
	ui.horizontalSlider->setSliderPosition(CurrentFrame);
}
void QTandOpenCV::setLayer(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	CurrentLayer = ui.setlayer->value();
	readCurrentFrame();
	showKeyFrameMark();
}
void QTandOpenCV::addLayer(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	layer b(FrameNum);
	pRoto->layers.push_back(b);
	CurrentLayer = pRoto->layers.size() - 1;
	ui.setlayer->setRange(0, pRoto->layers.size() - 1);
	ui.setlayer->setValue(CurrentLayer);
	ui.totalLayers->setText(QString::number(pRoto->layers.size()));

	//第0帧和FrameNum-1帧为关键帧
	pRoto->layers[CurrentLayer].KeyframeFlag[0] = true;
	pRoto->layers[CurrentLayer].KeyframeFlag[FrameNum - 1] = true;
	pRoto->layers[CurrentLayer].KeyframeIndex.push_back(0);
	pRoto->layers[CurrentLayer].KeyframeIndex.push_back(FrameNum - 1);

	ui.horizontalSlider->setValue(0);
	readCurrentFrame();
	showKeyFrameMark();
}
void QTandOpenCV::deleteLayer(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	vector<layer>::iterator it = pRoto->layers.begin() + CurrentLayer;
	pRoto->layers.erase(it);
	if (pRoto->layers.size() <= 0){
		layer b(FrameNum);
		pRoto->layers.push_back(b);
		CurrentLayer = 0;
		//第0帧和FrameNum-1帧为关键帧
		pRoto->layers[CurrentLayer].KeyframeFlag[0] = true;
		pRoto->layers[CurrentLayer].KeyframeFlag[FrameNum - 1] = true;
		pRoto->layers[CurrentLayer].KeyframeIndex.push_back(0);
		pRoto->layers[CurrentLayer].KeyframeIndex.push_back(FrameNum - 1);
	}
	else{
		if (CurrentLayer>0)
			CurrentLayer--;
	}
	ui.setlayer->setRange(0, pRoto->layers.size() - 1);
	ui.setlayer->setValue(CurrentLayer);
	ui.totalLayers->setText(QString::number(pRoto->layers.size()));
	ui.horizontalSlider->setValue(0);
	readCurrentFrame();
	showKeyFrameMark();
}
void QTandOpenCV::addKeyframe(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	if (ui.keyframe->isChecked()){
		pRoto->layers[CurrentLayer].KeyframeFlag[CurrentFrame] = true;
		if (pRoto->layers[CurrentLayer].KeyframeIndex.size() == 0 )  pRoto->layers[CurrentLayer].KeyframeIndex.push_back(CurrentFrame);
		else{
			vector<int>::iterator i = pRoto->layers[CurrentLayer].KeyframeIndex.begin();//迭代器
			vector<int>::iterator j = pRoto->layers[CurrentLayer].KeyframeIndex.end();//结尾，但是不包括
			while (i<j && *i<CurrentFrame){
				++i;
			}
			if (i == j)pRoto->layers[CurrentLayer].KeyframeIndex.push_back(CurrentFrame);
			else pRoto->layers[CurrentLayer].KeyframeIndex.insert(i, CurrentFrame);//在i前面插入一个元素
		}
		if (!pRoto->layers[CurrentLayer].GetResultDone){//如果当前帧还没有控制点了，就把最近关键帧的复制过来；如果得到结果后新增关键帧则已经有控制点了，所以不用复制。
			copyControlPoint();
			drawAllControlPoint();
			drawCurve();
			showFrame(frame);
		}
		
	}
	else{
		pRoto->layers[CurrentLayer].KeyframeFlag[CurrentFrame] = false;
		vector<int>::iterator i = pRoto->layers[CurrentLayer].KeyframeIndex.begin();
		while (*i<CurrentFrame){
			++i;
		}
		pRoto->layers[CurrentLayer].KeyframeIndex.erase(i);
	}
	showKeyFrameMark();
	if (pRoto->layers[CurrentLayer].GetResultDone)//更新相关帧的结果
		updateResult();
}

long long QTandOpenCV::imageTypeInSequence(string path, struct _finddata_t &fileinfo){//选择读取图片序列时候所需的类型
	long long hFile = -1;
	string p;
	hFile = (_findfirst(p.assign(path).append("\\*.PNG").c_str(), &fileinfo)); //默认搜索找到的第一种图片格式
	if (hFile == -1){
		hFile = (_findfirst(p.assign(path).append("\\*.JPG").c_str(), &fileinfo));
		if (hFile == -1){
			hFile = (_findfirst(p.assign(path).append("\\*.BMP").c_str(), &fileinfo));
		}
	}
	return hFile;
}

void QTandOpenCV::getFiles(string path, vector<string>& files) {
	//文件句柄  
	long long hFile = 0;//在win10中要改成long long
	//文件信息  
	struct _finddata_t fileinfo;   //大家可以去查看一下_finddata结构组成                            
	//以及_findfirst和_findnext的用法，了解后妈妈就再也不用担心我以后不会编了  
	string p;
	//if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1) {
	if ((hFile = imageTypeInSequence(path, fileinfo)) != -1) {//这样只遍历当前路径下的图片，而不会遍历子路径下的图片，且只读取多种类型中读到的第一种图片
		do {
			//如果是目录,迭代之  
			//如果不是,加入列表  
			if ((fileinfo.attrib & _A_SUBDIR)) {
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else {
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0); 
		_findclose(hFile);
	}
}

string QTandOpenCV::sequenceToVideo(string path){
	vector<string> files;
	getFiles(path, files);
	if (files.size() == 0) return "";
	Mat image = imread(files[0]);
	int width = image.cols;
	int height = image.rows;
	VideoWriter outputVideo;
	string VideoPath = path + "/sequenceToVideo.avi";
	outputVideo.open(VideoPath, CV_FOURCC('M', 'J', 'P', 'G'), 25.0, Size(width, height));
	if (!outputVideo.isOpened()) return "";

	outputVideo << image;
	for (size_t i = 1; i < files.size(); i++) {    //files.size()返回文件数量
		image = imread(files[i]);
		outputVideo << image;
	}
	return VideoPath;
}

//打开视频
void QTandOpenCV::openVideo()
{
	QMessageBox msgBox(QMessageBox::Warning, tr("input:"), tr("video or sequence?"), 0, this);
	msgBox.addButton(tr("video"), QMessageBox::AcceptRole);
	msgBox.addButton(tr("sequence"), QMessageBox::RejectRole);
	scale = 1.0;
	ui.scaleRatio->setText(QString::number(scale));
	ui.horizontalBar->setRange(0, 0);
	ui.verticalBar->setRange(0, 0);
	hTranslation = 0;
	vTranslation = 0;
	xPre = -1;
	yPre = -1;
	ControlPointFlag = Mat(height, width, CV_32SC1, Scalar(-1));//显示窗口的大小
	CurrentControlPoint = -1;

	if (msgBox.exec() == QMessageBox::AcceptRole){
		inputType = Video;

		QString file_full = QFileDialog::getOpenFileName(this, QString::fromLocal8Bit("请选择一个视频文件："), "D:/",tr("video files(*.avi *.mp4 *.rmvb *.wmv *.mkv )"));
		QFileInfo fileinfo = QFileInfo(file_full);
		QString fullfilepath = fileinfo.absoluteFilePath();//文件路径
		QString file_path = fileinfo.absolutePath();//绝对路径
		QString file_name = fileinfo.baseName();//文件名，不包括后缀
		videoPath = file_path.toStdString();//视频所在路径
		videoName = file_name.toStdString();//视频文件名
		if (!video.open(fullfilepath.toStdString())){//打开视频文件
			ui.label->setText("Cannot open the video!");
			return;
		}
		FrameNum = video.get(CV_CAP_PROP_FRAME_COUNT);
		Width = video.get(CV_CAP_PROP_FRAME_WIDTH);
		Height = video.get(CV_CAP_PROP_FRAME_HEIGHT);
	}
	else {
		inputType = Sequence;

		QString file_full = QFileDialog::getExistingDirectory(this, QString::fromLocal8Bit("请选择视频序列所在文件夹："), "D:/");
		if (file_full.isEmpty()) {
			ui.label->setText("Cannot load the sequence!");
			return;
		}
		videoPath = file_full.toStdString();
		videoName = "result";
		files.clear();//每次先清空一下
		getFiles((const char *)file_full.toLocal8Bit(), files);
		//需要一个函数来筛选图片文件
		
		if (files.size() == 0 ){//视频序列为空
			files.clear();//可能不是视频序列的文件夹，所以要清空一下读进来的文件列表
			ui.label->setText("no data in the sequence!");
			return;
		}
		//QString temp = QString::fromStdString(files[0]);
		//IplImage *ipl_image = cvLoadImage((const char *)temp.toLocal8Bit());
		//Mat ImageTemp(ipl_image, 0);//0不是复制，是公用data，各自的header
		Mat ImageTemp = imread(files[0]);//不能读有中文的路径
		FrameNum = files.size();
		Width = ImageTemp.cols;
		Height = ImageTemp.rows;
	}


	if (Height > height){
		ui.verticalBar->setRange(0,Height-height);
		ui.verticalBar->setSingleStep(50);
		ui.verticalBar->setPageStep(500);
	}
	if (Width > width){
		ui.horizontalBar->setRange(0, Width-width);
		ui.horizontalBar->setSingleStep(50);
		ui.horizontalBar->setPageStep(500);
	}
	

	scaledHeight = Height;
	scaledWidth = Width;
	CurrentFrame = 0;
	CurrentLayer = 0;
	ui.totalFrames->setText(QString::number(FrameNum));
	if (pRoto!=NULL){
		delete pRoto;
	}
	pRoto = new roto(FrameNum, Height, Width);
	ui.setlayer->setRange(0, pRoto->layers.size() - 1);
	ui.totalLayers->setText(QString::number(pRoto->layers.size()));

	

	//第0帧和FrameNum-1帧为关键帧
	pRoto->layers[CurrentLayer].KeyframeFlag[0] = true;
	pRoto->layers[CurrentLayer].KeyframeFlag[FrameNum - 1] = true;
	pRoto->layers[CurrentLayer].KeyframeIndex.push_back(0);
	pRoto->layers[CurrentLayer].KeyframeIndex.push_back(FrameNum - 1);


	ui.horizontalSlider->setMinimum(0);
	ui.horizontalSlider->setMaximum(FrameNum-1);
	ui.horizontalSlider->setValue(0);
	readCurrentFrame();
	showKeyFrameMark();
}

//显示一帧,输入的是已经resize之后的显示图像
void QTandOpenCV::showFrame(Mat input)
{
	Mat temp = input(Rect(hTranslation, vTranslation, min(scaledWidth, width), min(scaledHeight, height)));// (Size(width, height), CV_8UC3);
	//input(Rect(vTranslation, hTranslation, min(scaledWidth, width), min(scaledHeight, height))).convertTo(temp, pixelType, 1, 0);

	// 转换为QImage格式。QImage::Format_RGB888不同的摄像头用不同的格式。  
	// QImage(uchar *data, int width, int height, int bytesPerLine, Format format, QImageCleanupFunction cleanupFunction = 0, void *cleanupInfo = 0);
	QImage image = QImage(temp.data, temp.cols, temp.rows, temp.step, QImage::Format_RGB888).rgbSwapped();//指定每行的字节数，不然会出错
	ui.label->setPixmap(QPixmap::fromImage(image));  // 将图片显示到label上  
}

//读下一帧
void QTandOpenCV::readNextFrame()
{
	if (CurrentFrame >= FrameNum-1) return;//已经是最后一帧了

	if (inputType == Video){
		video >> frameOriginal;
		CurrentFrame = video.get(CV_CAP_PROP_POS_FRAMES);//出错
	}
	else if (inputType == Sequence){
		CurrentFrame++;
		frameOriginal = imread(files[CurrentFrame]);
	}
	
	ui.keyframe->setChecked(pRoto->layers[CurrentLayer].KeyframeFlag[CurrentFrame]);
	ui.currentFrame->setText(QString::number(CurrentFrame));
	ui.horizontalSlider->setSliderPosition(CurrentFrame);

	
	drawAllControlPointFlag();
	drawAllControlPoint();
	drawCurve();
	showFrame(frame);
}

//读当前帧
void QTandOpenCV::readCurrentFrame()
{
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size()==0) return;
	CurrentFrame = ui.horizontalSlider->value();
	ui.keyframe->setChecked(pRoto->layers[CurrentLayer].KeyframeFlag[CurrentFrame]);//死循环了?
	ui.currentFrame->setText(QString::number(CurrentFrame));


	if (inputType == Video){
		video.set(CV_CAP_PROP_POS_FRAMES, CurrentFrame);
		video >> frameOriginal;
	}
	else if (inputType == Sequence){
		frameOriginal = imread(files[CurrentFrame]);
	}

	
	drawAllControlPointFlag();//出问题?
	drawAllControlPoint();//已经有resize操作了
	drawCurve();
	showFrame(frame);
}



//析构函数
QTandOpenCV::~QTandOpenCV()
{
	if (inputType == Video){
		video.release();
	}
}





void QTandOpenCV::hMove(int value){
	hTranslation = value;
	drawAllControlPointFlag();
	drawAllControlPoint();
	drawCurve();
	showFrame(frame);
}


void QTandOpenCV::vMove(int value){
	vTranslation = value;
	drawAllControlPointFlag();
	drawAllControlPoint();
	drawCurve();
	showFrame(frame);
}

void QTandOpenCV::zoomInScale(){//记得更新拖动条的位置
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;

	if (scale >= 1.99) return;
	scale += 0.1;
	ui.scaleRatio->setText(QString::number(scale));
	scaledHeight = Height*scale;
	scaledWidth = Width*scale;
	if (scaledHeight > height){
		ui.verticalBar->setRange(0, scaledHeight - height);
	}
	else{
		ui.verticalBar->setRange(0, 0);
	}
	if (scaledWidth > width){
		ui.horizontalBar->setRange(0, scaledWidth - width);
	}
	else{
		ui.horizontalBar->setRange(0, 0);
	}

	drawAllControlPointFlag();
	drawAllControlPoint();
	drawCurve();
	showFrame(frame);
}

void QTandOpenCV::zoomOutScale(){
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;

	if (scale <= 0.51) return;
	scale -= 0.1;
	ui.scaleRatio->setText(QString::number(scale));
	scaledHeight = Height*scale;
	scaledWidth = Width*scale;


	//可能回超出范围，所以预先设置
	if (vTranslation>scaledHeight - height)
		vTranslation = max(0,scaledHeight - height);
	if (hTranslation > scaledWidth - width)
		hTranslation = max(0, scaledWidth - width);
	ui.verticalBar->setRange(0, max(0, scaledHeight - height));
	ui.horizontalBar->setRange(0, max(0, scaledWidth - width));

	
	drawAllControlPointFlag();
	drawAllControlPoint();
	drawCurve();
	showFrame(frame);
}


//--鼠标按下事件  
void QTandOpenCV::mousePressEvent(QMouseEvent *e)
{
	int x = e->x(), y = e->y();
	QString str = "mouse press: (" + QString::number(x) + ", " + QString::number(y) + ")";
	statusBar()->showMessage(str);
	
	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	if ((pRoto->layers[CurrentLayer].KeyframeFlag[CurrentFrame] || pRoto->layers[CurrentLayer].GetResultDone) && xyTest(x, y, fmin(scaledWidth, width), fmin(scaledHeight, height))){//是关键帧，在范围内，点击左键

		if (Qt::LeftButton == e->button()){//左键
			if (ControlPointFlag.at<int>(y, x) == -1){//还不是控制点
				if (CurrentFrame == 0 && !pRoto->layers[CurrentLayer].FirstKeyframeDone){//如果是第一帧，控制点还没完全添加，则添加控制点
					pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints.push_back(Point((x+hTranslation)/scale, (y+vTranslation)/scale));
					drawControlPointFlag(pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints.size() - 1, x, y);
					drawControlPoint(x, y);
					drawCurve();
					showFrame(frame);
				}
				else{//可能要平移所有控制点
					xPre = x;
					yPre = y;
				}
			}
			else{//已经是控制点了，则可能平移该控制点，或者最后一个控制点等于第一个控制点
				CurrentControlPoint = ControlPointFlag.at<int>(y, x);
				xPre = x;
				yPre = y;
			}
		}
		
		
	}

}

//---鼠标释放（松开）事件  
void QTandOpenCV::mouseReleaseEvent(QMouseEvent *e)
{
	int x = e->x(), y = e->y();
	int dx = x - xPre, dy = y - yPre;//平移量
	QString str = "mouse release: (" + QString::number(x) + ", " + QString::number(y) + ")";
	statusBar()->showMessage(str);

	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	int ControlPointSize = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints.size();
	if ((pRoto->layers[CurrentLayer].KeyframeFlag[CurrentFrame] || pRoto->layers[CurrentLayer].GetResultDone) && xyTest(x, y, fmin(scaledWidth, width), fmin(scaledHeight, height)) && (xPre != -1 && yPre != -1)){
		if (pRoto->layers[CurrentLayer].FirstKeyframeDone && CurrentControlPoint == -1  && (dx != 0 || dy != 0)){//第一个关键帧已经画完关键点了，点击的不是控制点，鼠标移动了，则平移所有控制点
			updateAllControlPoint(dx/scale,dy/scale);
			drawAllControlPoint();
			drawCurve();
			if (pRoto->layers[CurrentLayer].GetResultDone)//更新相关帧的结果
				updateResult();
		}
		else if (CurrentControlPoint != -1 && ControlPointFlag.at<int>(y, x) == -1){//点的时候是控制点，松开鼠标的时候不是控制点
			int ControlPointsSize = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints.size();
			if (CurrentControlPoint == 0 && ControlPointsSize>1 &&
				pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[0] == pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[ControlPointsSize-1]){//如果拖动的点事第一个控制点，且最后一个控制点也是同一点的话
				updateFirstControlPoint( dx/scale, dy/scale);//拖动第一个和最后一个点
				updateControlPoint(1, dx / scale, dy / scale);//同时拖动第二个点
				updateControlPoint(ControlPointSize-2, dx / scale, dy / scale);//同时拖动倒数第二个点
				drawAllControlPoint();
				drawCurve();
				if (pRoto->layers[CurrentLayer].GetResultDone)//更新相关帧的结果
					updateResult();
			}
			else{
				updateControlPoint(CurrentControlPoint, dx/scale, dy/scale);
				if (CurrentControlPoint%3==0){//bezier曲线的连接点，则平移相关控制点
					if (CurrentControlPoint - 1 > 0){
						updateControlPoint(CurrentControlPoint - 1, dx / scale, dy / scale);
					}
					if (CurrentControlPoint + 1 < ControlPointSize){
						updateControlPoint(CurrentControlPoint + 1, dx / scale, dy / scale);
					}
				}
				drawAllControlPoint();
				drawCurve();
				if (pRoto->layers[CurrentLayer].GetResultDone)//更新相关帧的结果
					updateResult();
			}
			
		}
		else if (CurrentFrame == 0 && ControlPointFlag.at<int>(y, x) == 0 && ControlPointFlag.at<int>(yPre, xPre) == 0 && ControlPointSize>0 && ControlPointSize % 3 == 0){//bezier曲线段的第四个点，且最后一个控制点等于第一个控制点
			Point pointTemp = pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints[0];
			pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints.push_back(pointTemp);
			drawCurve();

			pRoto->layers[CurrentLayer].FirstKeyframeDone = true;

			//复制第一个关键帧的控制点到其他关键帧
			/*for (int i = 1; i < pRoto->layers[CurrentLayer].KeyframeIndex.size(); ++i){//复制到其他关键帧
				pRoto->layers[CurrentLayer].frames[pRoto->layers[CurrentLayer].KeyframeIndex[i]].ControlPoints.assign(pRoto->layers[CurrentLayer].frames[0].ControlPoints.begin(),
					pRoto->layers[CurrentLayer].frames[0].ControlPoints.end());
			}*/
			//复制第一个关键帧的控制点到其他所有帧
			for (int i = 1; i < pRoto->FrameNum; ++i){
				pRoto->layers[CurrentLayer].frames[i].ControlPoints.assign(pRoto->layers[CurrentLayer].frames[0].ControlPoints.begin(),
					pRoto->layers[CurrentLayer].frames[0].ControlPoints.end());
			}
			pRoto->layers[CurrentLayer].GetResultDone = true;
		}
	}
	

	CurrentControlPoint = -1;
	xPre = -1;
	yPre = -1;
	showFrame(frame);
}

//--鼠标双击事件  
void QTandOpenCV::mouseDoubleClickEvent(QMouseEvent *e)
{
	//----Qt5解决中文乱码  
	QTextCodec *codec = QTextCodec::codecForName("GB18030");

	//----QMouseEvent类提供的x()和y()可获取鼠标相对窗口的位置  
	QString str = "(" + QString::number(e->x()) + ", " + QString::number(e->y()) + ")";
	statusBar()->showMessage(codec->toUnicode("鼠标双击:") + str, 1000);
}

//--鼠标移动事件  
void QTandOpenCV::mouseMoveEvent(QMouseEvent *e)
{
	//----QMouseEvent类提供的x()和y()可获取鼠标相对窗口的位置  
	QString str = "mouse move: (" + QString::number(e->x()) + ", " + QString::number(e->y()) + ")";
	statusBar()->showMessage(str);
}

Point QTandOpenCV::getCenter(vector<Point> ControlPoints){
	Point center(0, 0);
	for (int i = 0; i < ControlPoints.size();++i){
		center.x += ControlPoints[i].x;
		center.y += ControlPoints[i].y;
	}
	center.x = 1.0*center.x / ControlPoints.size() + 0.5;
	center.y = 1.0*center.y / ControlPoints.size() + 0.5;
	return center;
}

void QTandOpenCV::scaleControlPoints(vector<Point> & ControlPoints, float scale){
	Point center = getCenter(ControlPoints);
	for (int i = 0; i < ControlPoints.size(); ++i){
		ControlPoints[i].x = (ControlPoints[i].x - center.x)*scale + center.x;
		ControlPoints[i].y = (ControlPoints[i].y - center.y)*scale + center.y;
	}
}

//--鼠标滚动事件
void QTandOpenCV::wheelEvent(QWheelEvent * e){
	int x = e->x(), y = e->y();
	QString str = "mouse wheel move: (" + QString::number(x) + ", " + QString::number(y) + ")";
	statusBar()->showMessage(str);

	//如果视频还没开，则退出
	if (!video.isOpened() && files.size() == 0) return;
	if (pRoto->layers[CurrentLayer].FirstKeyframeDone && (pRoto->layers[CurrentLayer].KeyframeFlag[CurrentFrame] || pRoto->layers[CurrentLayer].GetResultDone) && xyTest(x, y, fmin(scaledWidth, width), fmin(scaledHeight, height))){
		// 当滚轮远离使用者时进行放大，当滚轮向使用者方向旋转时进行缩小
		if ((QApplication::keyboardModifiers() == Qt::ControlModifier)){//按下ctrl之后
			if (e->delta() > 0){
				scaleControlPoints(pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints, 1.05);
				drawAllControlPointFlag();
				drawAllControlPoint();
				drawCurve();
				showFrame(frame);
				if (pRoto->layers[CurrentLayer].GetResultDone)//更新相关帧的结果
					updateResult();
			}
			else if (e->delta() < 0){
				scaleControlPoints(pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints, 0.95);
				drawAllControlPointFlag();
				drawAllControlPoint();
				drawCurve();
				showFrame(frame);
				if (pRoto->layers[CurrentLayer].GetResultDone)//更新相关帧的结果
					updateResult();
			}
		}
		

	}
}

void QTandOpenCV::rotateControlPoints(vector<Point> & ControlPoints, float angle){
	Point center = getCenter(ControlPoints);
	for (int i = 0; i < ControlPoints.size(); ++i){
		int x_x0 = ControlPoints[i].x - center.x;
		int y_y0 = ControlPoints[i].y - center.y;
		float sinAngle = sin(angle);
		float cosAngle = cos(angle);
		ControlPoints[i].x = x_x0*cosAngle - y_y0*sinAngle + center.x;
		ControlPoints[i].y = x_x0*sinAngle + y_y0*cosAngle + center.y;
	}
}

void QTandOpenCV::keyPressEvent(QKeyEvent* e){
	if ((e->modifiers() == Qt::ControlModifier) && (e->key() == Qt::Key_A)){//逆时针旋转
		rotateControlPoints(pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints, 2 * PI * 0.99);;
		drawAllControlPointFlag();
		drawAllControlPoint();
		drawCurve();
		showFrame(frame);
		if (pRoto->layers[CurrentLayer].GetResultDone)//更新相关帧的结果
			updateResult();
	}
	else if ((e->modifiers() == Qt::ControlModifier) && (e->key() == Qt::Key_C)){//顺时针旋转
		rotateControlPoints(pRoto->layers[CurrentLayer].frames[CurrentFrame].ControlPoints, 2 * PI * 0.01);
		drawAllControlPointFlag();
		drawAllControlPoint();
		drawCurve();
		showFrame(frame);
		if (pRoto->layers[CurrentLayer].GetResultDone)//更新相关帧的结果
			updateResult();
	}
}
