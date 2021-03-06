SOURCES += loghandler.cpp \
    $$PWD/fid.cpp \
    $$PWD/ftworker.cpp \
    $$PWD/experiment.cpp \
    $$PWD/ftmwconfig.cpp \
    $$PWD/chirpconfig.cpp \
    $$PWD/pulsegenconfig.cpp \
    $$PWD/flowconfig.cpp \
    $$PWD/lifconfig.cpp \
    $$PWD/liftrace.cpp \
    $$PWD/datastructs.cpp \
    $$PWD/analysis.cpp \
    $$PWD/ioboardconfig.cpp \
    $$PWD/snapworker.cpp \
    $$PWD/peakfinder.cpp

HEADERS += loghandler.h \
    $$PWD/fid.h \
    $$PWD/ftworker.h \
    $$PWD/experiment.h \
    $$PWD/ftmwconfig.h \
    $$PWD/gpuaverager.h \
    $$PWD/chirpconfig.h \
    $$PWD/pulsegenconfig.h \
    $$PWD/flowconfig.h \
    $$PWD/lifconfig.h \
    $$PWD/liftrace.h \
    $$PWD/datastructs.h \
    $$PWD/analysis.h \
    $$PWD/ioboardconfig.h \
    $$PWD/snapworker.h \
    $$PWD/peakfinder.h

gpu-cuda {
    OTHER_FILES += $$PWD/gpuaverager.cu
}
