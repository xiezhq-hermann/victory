FC          = mpifort
FC_FLAGS    = -traceback -Mbounds -ta=tesla:managed -Minfo=accel
DIR         = /home/richard/OpenACC_Hackathon/victory
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
