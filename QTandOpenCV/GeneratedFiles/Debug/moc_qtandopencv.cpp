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
    QByteArrayData data[12];
    char stringdata[127];
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
QT_MOC_LITERAL(3, 23, 9),
QT_MOC_LITERAL(4, 33, 16),
QT_MOC_LITERAL(5, 50, 9),
QT_MOC_LITERAL(6, 60, 10),
QT_MOC_LITERAL(7, 71, 11),
QT_MOC_LITERAL(8, 83, 8),
QT_MOC_LITERAL(9, 92, 8),
QT_MOC_LITERAL(10, 101, 12),
QT_MOC_LITERAL(11, 114, 12)
    },
    "QTandOpenCV\0openVideo\0\0readFrame\0"
    "readCurrentFrame\0stopVideo\0startVideo\0"
    "addKeyframe\0addLayer\0setLayer\0"
    "lastKeyframe\0nextKeyframe"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_QTandOpenCV[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      10,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   64,    2, 0x08 /* Private */,
       3,    0,   65,    2, 0x08 /* Private */,
       4,    0,   66,    2, 0x08 /* Private */,
       5,    0,   67,    2, 0x08 /* Private */,
       6,    0,   68,    2, 0x08 /* Private */,
       7,    0,   69,    2, 0x08 /* Private */,
       8,    0,   70,    2, 0x08 /* Private */,
       9,    0,   71,    2, 0x08 /* Private */,
      10,    0,   72,    2, 0x08 /* Private */,
      11,    0,   73,    2, 0x08 /* Private */,

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

       0        // eod
};

void QTandOpenCV::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        QTandOpenCV *_t = static_cast<QTandOpenCV *>(_o);
        switch (_id) {
        case 0: _t->openVideo(); break;
        case 1: _t->readFrame(); break;
        case 2: _t->readCurrentFrame(); break;
        case 3: _t->stopVideo(); break;
        case 4: _t->startVideo(); break;
        case 5: _t->addKeyframe(); break;
        case 6: _t->addLayer(); break;
        case 7: _t->setLayer(); break;
        case 8: _t->lastKeyframe(); break;
        case 9: _t->nextKeyframe(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
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
        if (_id < 10)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 10;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 10)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 10;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
