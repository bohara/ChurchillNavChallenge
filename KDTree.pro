#-------------------------------------------------
#
# Project created by QtCreator 2015-03-02T22:21:43
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = KDTree
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

CONFIG += c++11

SOURCES += main.cpp \
    kdtree.cpp

HEADERS += \
    kdtree.h \
    point.h \
    boundedpqueue.h
