TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    integrator.cpp

HEADERS += \

CXX_FLAGS += -O3 -std=c++0x

INCLUDEPATH += /home/goranbs/goran/CompPhys/programs/cppLibrary/
