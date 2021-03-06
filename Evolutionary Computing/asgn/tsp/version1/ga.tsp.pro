######################################################################
# Automatically generated by qmake (3.0) Fri Mar 13 20:18:25 2015
######################################################################

TEMPLATE             = app
TARGET               = ga.tsp
CONFIG              += c++11
DEFINES             += __DEBUG__ __QCREATOR__
INCLUDEPATH         += . inc vendors/jsoncons/src
# General libs
unix:LIBS           += -lpthread -lboost_program_options -lboost_regex -lboost_system -lboost_filesystem -lboost_serialization

# Input
HEADERS += inc/chromosome.hpp \
    inc/stdafx.hpp \
    inc/population.hpp \
    inc/gaConfig.hpp \
    inc/utils.hpp \
    inc/data.hpp \
    main.helper.hpp \
    inc/Timer.hpp
SOURCES += main.cpp \
    src/population.cpp \

OTHER_FILES += \
    config.json \
    dependency.json
