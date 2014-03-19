TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    integrator.cpp \
    unitconverter.cpp \
    cell.cpp \
    initialstate.cpp \
    atom.cpp \
    thermostat.cpp

#LIBS += -larmadillo -lblas -llapack

HEADERS += \
    unitconverter.h \
    cell.h \
    initialstate.h \
    atom.h

release {
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

INCLUDEPATH += /home/goranbs/goran/CompPhys/programs/cppLibrary/
