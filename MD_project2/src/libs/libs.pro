include(../../defaults.pri)
TEMPLATE = lib
TARGET= myapp

SOURCES = atom.cpp \
      unitconverter.cpp \
      initialstate.cpp \
    generatenanoporoussystem.cpp \
    readinitialstatefile.cpp

HEADERS = atom.h \
      unitconverter.h \
      initialstate.h \
    generatenanoporoussystem.h \
    readinitialstatefile.h

LIBS += -larmadillo -lblas -llapack
