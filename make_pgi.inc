FC          = mpifort
FC_FLAGS    = -ta=tesla:cc60 -Minfo=accel -Mcuda -Mcudalib=cufft -Mcudalib=cublas
DIR         = /home/gavin/Documents/github/victory
BIN         = $(DIR)/bin
FLIB        = -lblas -llapack -lpthread
MODPATH     = -I$(DIR)/src/SRC_MOD

LLIB = $(DIR)/lib/lib.a $(DIR)/lib/libdfftpack.a

LIB = $(LLIB) $(FLIB)

#
#  Archiver and its flags
#

ARCH   = ar
ARFLAG = cr
RANLIB = ranlib

# set flags
FLAGS = $(FC_FLAGS)
