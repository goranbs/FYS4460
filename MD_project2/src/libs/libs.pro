include(../../defaults.pri)
TEMPLATE = lib
TARGET= myapp

SOURCES = atom.cpp \
      unitconverter.cpp \
      initialstate.cpp

HEADERS = atom.h \
      unitconverter.h \
      initialstate.h

LIBS += -larmadillo -lblas -llapack
