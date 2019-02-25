QT += core
QT -= gui

TARGET = 1DElasticProject
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    statevector.cpp \
    equationofstate.cpp \
    forcesolver.cpp \
    vectoralgebra.cpp \
    solvers.cpp \
    matrixalgebra.cpp \
    tensoralgebra.cpp \
    elasticequationofstate.cpp \
    hyperelasticvariables.cpp \
    elasticstatevector.cpp \
    elasticsolvers.cpp \
    elasticforcesolver.cpp \
    slopelimiters.cpp \
    slicsolver.cpp \
    elasticslopelimiters.cpp \
    elasticslicsolver.cpp \
    tests.cpp \
    unittests.cpp

HEADERS += \
    statevector.h \
    equationofstate.h \
    forcesolver.h \
    vectoralgebra.h \
    solvers.h \
    matrixalgebra.h \
    tensoralgebra.h \
    elasticequationofstate.h \
    hyperelasticvariables.h \
    elasticstatevector.h \
    elasticsolvers.h \
    elasticforcesolver.h \
    slopelimiters.h \
    slicsolver.h \
    elasticslopelimiters.h \
    elasticslicsolver.h \
    tests.h \
    unittests.h

