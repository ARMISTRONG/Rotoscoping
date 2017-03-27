/****************************************************************************
** Meta object code from reading C++ file 'qtandopencv.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.3.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../qtandopencv.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'qtandopencv.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.3.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_QTandOpenCV_t {
    QByteArrayData data[23];
    char stringdata[277];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_QTandOpenCV_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_QTandOpenCV_t qt_meta_stringdata_QTandOpenCV = {
    {
QT_MOC_LITERAL(0, 0, 11),
QT_MOC_LITERAL(1, 12, 9),
QT_MOC_LITERAL(2, 22, 0),
QT_MOC_LITERAL(3, 23, 14),
QT_MOC_LITERAL(4, 38, 13),
QT_MOC_LITERAL(5, 52, 16),
QT_MOC_LITERAL(6, 69, 11),
QT_MOC_LITERAL(7, 81, 8),
QT_MOC_LITERAL(8, 90, 11),
QT_MOC_LITERAL(9, 102, 8),
QT_MOC_LITERAL(10, 111, 12),
QT_MOC_LITERAL(11, 124, 12),
QT_MOC_LITERAL(12, 137, 10),
QT_MOC_LITERAL(13, 148, 9),
QT_MOC_LITERAL(14, 158, 12),
QT_MOC_LITERAL(15, 171, 15),
QT_MOC_LITERAL(16, 187, 17),
QT_MOC_LITERAL(17, 205, 16),
QT_MOC_LITERAL(18, 222, 17),
QT_MOC_LITERAL(19, 240, 5),
QT_MOC_LITERAL(20, 246, 5),
QT_MOC_LITERAL(21, 252, 11),
QT_MOC_LITERAL(22, 264, 12)
    },
    "QTandOpenCV\0openVideo\0\0readResultText\0"
    "readNextFrame\0readCurrentFrame\0"
    "addKeyframe\0addLayer\0deleteLayer\0"
    "setLayer\0lastKeyframe\0nextKeyframe\0"
    "saveResult\0getResult\0setLinearFit\0"
    "setParabolicFit\0setCubicSplineFit\0"
    "setPolynomialFit\0setPolynomialFitN\0"
    "hMove\0vMove\0zoomInScale\0zoomOutScale"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_QTandOpenCV[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      21,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  119,    2, 0x08 /* Private */,
       3,    0,  120,    2, 0x08 /* Private */,
       4,    0,  121,    2, 0x08 /* Private */,
       5,    0,  122,    2, 0x08 /* Private */,
       6,    0,  123,    2, 0x08 /* Private */,
       7,    0,  124,    2, 0x08 /* Private */,
       8,    0,  125,    2, 0x08 /* Private */,
       9,    0,  126,    2, 0x08 /* Private */,
      10,    0,  127,    2, 0x08 /* Private */,
      11,    0,  128,    2, 0x08 /* Private */,
      12,    0,  129,    2, 0x08 /* Private */,
      13,    0,  130,    2, 0x08 /* Private */,
      14,    0,  131,    2, 0x08 /* Private */,
      15,    0,  132,    2, 0x08 /* Private */,
      16,    0,  133,    2, 0x08 /* Private */,
      17,    0,  134,    2, 0x08 /* Private */,
      18,    0,  135,    2, 0x08 /* Private */,
      19,    1,  136,    2, 0x08 /* Private */,
      20,    1,  139,    2, 0x08 /* Private */,
      21,    0,  142,    2, 0x08 /* Private */,
      22,    0,  143,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void QTandOpenCV::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        QTandOpenCV *_t = static_cast<QTandOpenCV *>(_o);
        switch (_id) {
        case 0: _t->openVideo(); break;
        case 1: _t->readResultText(); break;
        case 2: _t->readNextFrame(); break;
        case 3: _t->readCurrentFrame(); break;
        case 4: _t->addKeyframe(); break;
        case 5: _t->addLayer(); break;
        case 6: _t->deleteLayer(); break;
        case 7: _t->setLayer(); break;
        case 8: _t->lastKeyframe(); break;
        case 9: _t->nextKeyframe(); break;
        case 10: _t->saveResult(); break;
        case 11: _t->getResult(); break;
        case 12: _t->setLinearFit(); break;
        case 13: _t->setParabolicFit(); break;
        case 14: _t->setCubicSplineFit(); break;
        case 15: _t->setPolynomialFit(); break;
        case 16: _t->setPolynomialFitN(); break;
        case 17: _t->hMove((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 18: _t->vMove((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 19: _t->zoomInScale(); break;
        case 20: _t->zoomOutScale(); break;
        default: ;
        }
    }
}

const QMetaObject QTandOpenCV::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_QTandOpenCV.data,
      qt_meta_data_QTandOpenCV,  qt_static_metacall, 0, 0}
};


const QMetaObject *QTandOpenCV::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *QTandOpenCV::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_QTandOpenCV.stringdata))
        return static_cast<void*>(const_cast< QTandOpenCV*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int QTandOpenCV::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 21)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 21;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 21)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 21;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
