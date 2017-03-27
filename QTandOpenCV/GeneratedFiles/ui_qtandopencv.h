/********************************************************************************
** Form generated from reading UI file 'qtandopencv.ui'
**
** Created by: Qt User Interface Compiler version 5.3.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QTANDOPENCV_H
#define UI_QTANDOPENCV_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QScrollBar>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_QTandOpenCVClass
{
public:
    QWidget *centralWidget;
    QLabel *label;
    QSlider *horizontalSlider;
    QGroupBox *interpolation;
    QRadioButton *linearfit;
    QRadioButton *polynomialfit;
    QSpinBox *npolynomialfit;
    QRadioButton *parabolicfit;
    QRadioButton *cubicsplinefit;
    QGroupBox *result;
    QPushButton *getresult;
    QPushButton *saveresult;
    QScrollBar *horizontalBar;
    QScrollBar *verticalBar;
    QGroupBox *layer;
    QLabel *totalLayers;
    QLabel *label_8;
    QSpinBox *setlayer;
    QLabel *label_7;
    QPushButton *addlayer;
    QPushButton *deletelayer;
    QGroupBox *scaling;
    QPushButton *zoomOut;
    QPushButton *zoomIn;
    QLabel *scaleRatio;
    QLabel *label_9;
    QGroupBox *frame;
    QPushButton *nextkeyframe;
    QPushButton *lastkeyframe;
    QLabel *totalFrames;
    QLabel *currentFrame;
    QLabel *label_5;
    QLabel *label_4;
    QRadioButton *keyframe;
    QGroupBox *load;
    QPushButton *open;
    QPushButton *loadresult;
    QLabel *keyframemark;
    QMenuBar *menuBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *QTandOpenCVClass)
    {
        if (QTandOpenCVClass->objectName().isEmpty())
            QTandOpenCVClass->setObjectName(QStringLiteral("QTandOpenCVClass"));
        QTandOpenCVClass->setEnabled(true);
        QTandOpenCVClass->resize(1296, 960);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(QTandOpenCVClass->sizePolicy().hasHeightForWidth());
        QTandOpenCVClass->setSizePolicy(sizePolicy);
        QTandOpenCVClass->setMinimumSize(QSize(1296, 960));
        QTandOpenCVClass->setMaximumSize(QSize(1296, 960));
        QPalette palette;
        QBrush brush(QColor(255, 255, 255, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Base, brush);
        QBrush brush1(QColor(100, 100, 100, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Window, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Base, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Window, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Window, brush1);
        QTandOpenCVClass->setPalette(palette);
        QTandOpenCVClass->setStyleSheet(QStringLiteral(""));
        centralWidget = new QWidget(QTandOpenCVClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        sizePolicy.setHeightForWidth(centralWidget->sizePolicy().hasHeightForWidth());
        centralWidget->setSizePolicy(sizePolicy);
        label = new QLabel(centralWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(0, 0, 1280, 720));
        sizePolicy.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy);
        label->setMinimumSize(QSize(1280, 720));
        label->setMaximumSize(QSize(1280, 720));
        label->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        horizontalSlider = new QSlider(centralWidget);
        horizontalSlider->setObjectName(QStringLiteral("horizontalSlider"));
        horizontalSlider->setGeometry(QRect(0, 736, 1296, 22));
        horizontalSlider->setOrientation(Qt::Horizontal);
        interpolation = new QGroupBox(centralWidget);
        interpolation->setObjectName(QStringLiteral("interpolation"));
        interpolation->setGeometry(QRect(880, 780, 161, 151));
        linearfit = new QRadioButton(interpolation);
        linearfit->setObjectName(QStringLiteral("linearfit"));
        linearfit->setGeometry(QRect(10, 25, 111, 16));
        sizePolicy.setHeightForWidth(linearfit->sizePolicy().hasHeightForWidth());
        linearfit->setSizePolicy(sizePolicy);
        linearfit->setContextMenuPolicy(Qt::DefaultContextMenu);
        linearfit->setChecked(true);
        polynomialfit = new QRadioButton(interpolation);
        polynomialfit->setObjectName(QStringLiteral("polynomialfit"));
        polynomialfit->setGeometry(QRect(10, 115, 121, 20));
        sizePolicy.setHeightForWidth(polynomialfit->sizePolicy().hasHeightForWidth());
        polynomialfit->setSizePolicy(sizePolicy);
        polynomialfit->setCheckable(false);
        npolynomialfit = new QSpinBox(interpolation);
        npolynomialfit->setObjectName(QStringLiteral("npolynomialfit"));
        npolynomialfit->setGeometry(QRect(100, 115, 42, 22));
        parabolicfit = new QRadioButton(interpolation);
        parabolicfit->setObjectName(QStringLiteral("parabolicfit"));
        parabolicfit->setGeometry(QRect(10, 55, 111, 16));
        parabolicfit->setCheckable(false);
        cubicsplinefit = new QRadioButton(interpolation);
        cubicsplinefit->setObjectName(QStringLiteral("cubicsplinefit"));
        cubicsplinefit->setGeometry(QRect(10, 85, 111, 16));
        result = new QGroupBox(centralWidget);
        result->setObjectName(QStringLiteral("result"));
        result->setGeometry(QRect(1070, 780, 181, 61));
        getresult = new QPushButton(result);
        getresult->setObjectName(QStringLiteral("getresult"));
        getresult->setGeometry(QRect(10, 20, 71, 23));
        saveresult = new QPushButton(result);
        saveresult->setObjectName(QStringLiteral("saveresult"));
        saveresult->setGeometry(QRect(90, 20, 75, 23));
        horizontalBar = new QScrollBar(centralWidget);
        horizontalBar->setObjectName(QStringLiteral("horizontalBar"));
        horizontalBar->setGeometry(QRect(0, 720, 1280, 16));
        horizontalBar->setMaximum(0);
        horizontalBar->setOrientation(Qt::Horizontal);
        verticalBar = new QScrollBar(centralWidget);
        verticalBar->setObjectName(QStringLiteral("verticalBar"));
        verticalBar->setGeometry(QRect(1280, 0, 16, 736));
        verticalBar->setMaximum(0);
        verticalBar->setOrientation(Qt::Vertical);
        layer = new QGroupBox(centralWidget);
        layer->setObjectName(QStringLiteral("layer"));
        layer->setGeometry(QRect(680, 780, 171, 131));
        totalLayers = new QLabel(layer);
        totalLayers->setObjectName(QStringLiteral("totalLayers"));
        totalLayers->setGeometry(QRect(110, 60, 51, 16));
        label_8 = new QLabel(layer);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(5, 25, 101, 20));
        sizePolicy.setHeightForWidth(label_8->sizePolicy().hasHeightForWidth());
        label_8->setSizePolicy(sizePolicy);
        label_8->setAlignment(Qt::AlignCenter);
        setlayer = new QSpinBox(layer);
        setlayer->setObjectName(QStringLiteral("setlayer"));
        setlayer->setGeometry(QRect(110, 25, 41, 22));
        label_7 = new QLabel(layer);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setGeometry(QRect(5, 55, 95, 20));
        label_7->setAlignment(Qt::AlignCenter);
        addlayer = new QPushButton(layer);
        addlayer->setObjectName(QStringLiteral("addlayer"));
        addlayer->setGeometry(QRect(15, 90, 61, 23));
        deletelayer = new QPushButton(layer);
        deletelayer->setObjectName(QStringLiteral("deletelayer"));
        deletelayer->setGeometry(QRect(80, 90, 81, 23));
        scaling = new QGroupBox(centralWidget);
        scaling->setObjectName(QStringLiteral("scaling"));
        scaling->setGeometry(QRect(240, 780, 161, 91));
        zoomOut = new QPushButton(scaling);
        zoomOut->setObjectName(QStringLiteral("zoomOut"));
        zoomOut->setGeometry(QRect(10, 50, 61, 23));
        zoomIn = new QPushButton(scaling);
        zoomIn->setObjectName(QStringLiteral("zoomIn"));
        zoomIn->setGeometry(QRect(80, 50, 61, 23));
        scaleRatio = new QLabel(scaling);
        scaleRatio->setObjectName(QStringLiteral("scaleRatio"));
        scaleRatio->setGeometry(QRect(90, 25, 31, 16));
        label_9 = new QLabel(scaling);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setGeometry(QRect(0, 20, 95, 20));
        label_9->setAlignment(Qt::AlignCenter);
        frame = new QGroupBox(centralWidget);
        frame->setObjectName(QStringLiteral("frame"));
        frame->setGeometry(QRect(430, 780, 221, 151));
        nextkeyframe = new QPushButton(frame);
        nextkeyframe->setObjectName(QStringLiteral("nextkeyframe"));
        nextkeyframe->setGeometry(QRect(110, 110, 91, 23));
        lastkeyframe = new QPushButton(frame);
        lastkeyframe->setObjectName(QStringLiteral("lastkeyframe"));
        lastkeyframe->setGeometry(QRect(10, 110, 91, 23));
        totalFrames = new QLabel(frame);
        totalFrames->setObjectName(QStringLiteral("totalFrames"));
        totalFrames->setGeometry(QRect(100, 80, 51, 16));
        currentFrame = new QLabel(frame);
        currentFrame->setObjectName(QStringLiteral("currentFrame"));
        currentFrame->setGeometry(QRect(100, 55, 51, 16));
        currentFrame->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);
        label_5 = new QLabel(frame);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(10, 75, 81, 20));
        label_5->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);
        label_4 = new QLabel(frame);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(10, 50, 91, 20));
        label_4->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);
        keyframe = new QRadioButton(frame);
        keyframe->setObjectName(QStringLiteral("keyframe"));
        keyframe->setGeometry(QRect(10, 25, 71, 16));
        load = new QGroupBox(centralWidget);
        load->setObjectName(QStringLiteral("load"));
        load->setGeometry(QRect(30, 780, 181, 61));
        open = new QPushButton(load);
        open->setObjectName(QStringLiteral("open"));
        open->setGeometry(QRect(10, 20, 71, 23));
        loadresult = new QPushButton(load);
        loadresult->setObjectName(QStringLiteral("loadresult"));
        loadresult->setGeometry(QRect(90, 20, 75, 23));
        keyframemark = new QLabel(centralWidget);
        keyframemark->setObjectName(QStringLiteral("keyframemark"));
        keyframemark->setGeometry(QRect(0, 745, 1296, 21));
        keyframemark->setAlignment(Qt::AlignCenter);
        QTandOpenCVClass->setCentralWidget(centralWidget);
        keyframemark->raise();
        label->raise();
        interpolation->raise();
        result->raise();
        horizontalBar->raise();
        verticalBar->raise();
        layer->raise();
        scaling->raise();
        frame->raise();
        currentFrame->raise();
        load->raise();
        horizontalSlider->raise();
        menuBar = new QMenuBar(QTandOpenCVClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1296, 23));
        QTandOpenCVClass->setMenuBar(menuBar);
        statusBar = new QStatusBar(QTandOpenCVClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        QTandOpenCVClass->setStatusBar(statusBar);

        retranslateUi(QTandOpenCVClass);

        QMetaObject::connectSlotsByName(QTandOpenCVClass);
    } // setupUi

    void retranslateUi(QMainWindow *QTandOpenCVClass)
    {
        QTandOpenCVClass->setWindowTitle(QApplication::translate("QTandOpenCVClass", "QTandOpenCV", 0));
        label->setText(QString());
        interpolation->setTitle(QApplication::translate("QTandOpenCVClass", "interpolation", 0));
        linearfit->setText(QApplication::translate("QTandOpenCVClass", "linear", 0));
        polynomialfit->setText(QApplication::translate("QTandOpenCVClass", "polynomial:", 0));
        parabolicfit->setText(QApplication::translate("QTandOpenCVClass", "parabolic", 0));
        cubicsplinefit->setText(QApplication::translate("QTandOpenCVClass", "cubicspline", 0));
        result->setTitle(QApplication::translate("QTandOpenCVClass", "result", 0));
        getresult->setText(QApplication::translate("QTandOpenCVClass", "get result", 0));
        saveresult->setText(QApplication::translate("QTandOpenCVClass", "save result", 0));
        layer->setTitle(QApplication::translate("QTandOpenCVClass", "layer", 0));
        totalLayers->setText(QApplication::translate("QTandOpenCVClass", "0", 0));
        label_8->setText(QApplication::translate("QTandOpenCVClass", "current layer:", 0));
        label_7->setText(QApplication::translate("QTandOpenCVClass", "total layers:", 0));
        addlayer->setText(QApplication::translate("QTandOpenCVClass", "add layer", 0));
        deletelayer->setText(QApplication::translate("QTandOpenCVClass", "delete layer", 0));
        scaling->setTitle(QApplication::translate("QTandOpenCVClass", "scaling", 0));
        zoomOut->setText(QApplication::translate("QTandOpenCVClass", "zoom out", 0));
        zoomIn->setText(QApplication::translate("QTandOpenCVClass", "zoom in", 0));
        scaleRatio->setText(QApplication::translate("QTandOpenCVClass", "0", 0));
        label_9->setText(QApplication::translate("QTandOpenCVClass", "scale ratio:", 0));
        frame->setTitle(QApplication::translate("QTandOpenCVClass", "frame", 0));
        nextkeyframe->setText(QApplication::translate("QTandOpenCVClass", "next keyframe", 0));
        lastkeyframe->setText(QApplication::translate("QTandOpenCVClass", "last keyframe", 0));
        totalFrames->setText(QApplication::translate("QTandOpenCVClass", "0", 0));
        currentFrame->setText(QApplication::translate("QTandOpenCVClass", "0", 0));
        label_5->setText(QApplication::translate("QTandOpenCVClass", "total frames:", 0));
        label_4->setText(QApplication::translate("QTandOpenCVClass", "current frame:", 0));
        keyframe->setText(QApplication::translate("QTandOpenCVClass", "keyframe", 0));
        load->setTitle(QApplication::translate("QTandOpenCVClass", "load", 0));
        open->setText(QApplication::translate("QTandOpenCVClass", "load video", 0));
        loadresult->setText(QApplication::translate("QTandOpenCVClass", "load result", 0));
        keyframemark->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class QTandOpenCVClass: public Ui_QTandOpenCVClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QTANDOPENCV_H
