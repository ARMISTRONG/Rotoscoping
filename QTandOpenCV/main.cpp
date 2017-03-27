#include "qtandopencv.h"
#include <QtWidgets/QApplication>
#include <QLabel>
#include <QPushButton>
#include <QTextCodec>
#include <head.h>


int main(int argc, char *argv[])
{
	QTextCodec::setCodecForLocale(QTextCodec::codecForName("GBK"));


	QApplication a(argc, argv);
	QTandOpenCV w;
	w.show();
	return a.exec();
}