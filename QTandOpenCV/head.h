#ifndef Head_H
#define Head_H

#include <iostream>
#include <fstream>//使用文件流进行操作，主要使用的类是ifstream, ofstream
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <direct.h>//获取当前路径
#include <time.h>
#include <iomanip>//设置输出精度
//#include <windows.h>//与opencv冲突
#include <thread>//多线程
#include <mutex>//同步互斥
#include <limits>



#include "cv.h"
#include "cxcore.h"
#include "cvaux.h"
#include "highgui.h"
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

//和sift特征提取有关的库
#include <opencv2/features2d/features2d.hpp>
#include<opencv2/nonfree/nonfree.hpp>
#include<opencv2/legacy/legacy.hpp>

//读取目录下所有图片
#include <iostream>  
#include <string>  
#include <vector>  
#include <io.h>  
#include <string.h>  

using namespace std;
using namespace cv;


#endif