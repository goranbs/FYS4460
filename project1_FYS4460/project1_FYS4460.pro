TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    integrator.cpp \
    unitconverter.cpp \
    cell.cpp

HEADERS += \
    unitconverter.h \
    cell.h

QMAKE_CXXFLAGS += -O3 -std=c++0x

INCLUDEPATH += /home/goranbs/goran/CompPhys/programs/cppLibrary/
