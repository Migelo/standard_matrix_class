QT -= core
QT -= gui

QMAKE_CXXFLAGS -= -std=gnu++98
QMAKE_CXXFLAGS -= -std=gnu++0x
QMAKE_CXXFLAGS -= --std=c++11
QMAKE_CXXFLAGS_RELEASE -= -std=gnu++98
QMAKE_CXXFLAGS_RELEASE -= -std=gnu++0x
QMAKE_CXXFLAGS_RELEASE -= --std=c++11

QMAKE_CXXFLAGS += -std=c++14

TARGET = aacp_template
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

CONFIG += c++11

LIBS += -L/usr/lib -llapack
LIBS += -L/usr/lib -lblas

INCLUDEPATH += /opt/local/include

SOURCES += main.cpp \
    la_objects.cpp

HEADERS += \
    la_wrapper.h \
    la_base_obj.h \
    la_objects.h \
    la_operations.h

