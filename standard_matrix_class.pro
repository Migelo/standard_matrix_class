QT -= core
QT -= gui

QMAKE_CXXFLAGS -= -std=gnu++98
QMAKE_CXXFLAGS -= -std=gnu++0x
QMAKE_CXXFLAGS -= --std=c++11
QMAKE_CXXFLAGS_RELEASE -= -std=gnu++98
QMAKE_CXXFLAGS_RELEASE -= -std=gnu++0x
QMAKE_CXXFLAGS_RELEASE -= --std=c++11

TARGET = aacp_template
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

CONFIG += c++11

LIBS += -L/usr/lib -llapack
LIBS += -L/usr/lib -lblas
#LIBS += -framework Accelerate
INCLUDEPATH += /usr/include
#INCLUDEPATH += /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers/

SOURCES += main.cpp \
    la_objects.cpp

HEADERS += \
    la_wrapper.h \
    la_base_obj.h \
    la_objects.h \
    la_operations.h \
    lapack_header.h \
    my_operations.h \
    my_ops_wrapper.h

