TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11

INCLUDEPATH += $$system(${ROOTSYS}/bin/root-config --incdir)
INCLUDEPATH += $$system(${ROOTSYS}/bin/root-config --libs)
DEPENDPATH  += $$system(${ROOTSYS}/bin/root-config --libs)
LIBS        += $$system(${ROOTSYS}/bin/root-config --libs)
LIBS        += $$system(${ROOTSYS}/bin/root-config --glibs)

TARGET = ../Scripts/entanglement

SOURCES += main.cpp \
    analysis.cpp \
    entangled.cpp\
    video.cpp \
    dual.cpp

HEADERS += \
    analysis.h \
    entangled.h\
    video.h \
    dual.h
