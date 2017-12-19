# Hey emacs, this is a -*- makefile -*-

# Copyright (c) 2015 Open Source Medical Software Corporation,
#                    University of California, San Diego.
# 
# All rights reserved. 
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, merge, publish, 
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
# 
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

TOP=../../../../../BuildWithMake

include $(TOP)/include.mk

CCFLAGS = $(GLOBAL_CCFLAGS) $(MPI_INCDIR) -I .

HDRS	= defs.h macros.h metis.h \
          proto.h rename.h struct.h

CSRCS	= comm.c util.c debug.c setup.c grsetup.c timer.c \
       node_refine.c initmsection.c \
       xyzpart.c pspases.c frename.c \
       iintsort.c iidxsort.c ikeysort.c ikeyvalsort.c \
       kmetis.c gkmetis.c ometis.c \
       initpart.c match.c \
       kwayfm.c kwayrefine.c kwaybalance.c \
       remap.c stat.c fpqueue.c \
       ametis.c rmetis.c initbalance.c \
       mdiffusion.c diffutil.c wave.c \
       csrmatch.c redomylink.c balancemylink.c \
       selectq.c akwayfm.c serial.c move.c \
       mmetis.c mesh.c memory.c weird.c backcompat.c

TARGET_LIB_NAME = simvascular_thirdparty_parmetis

include $(TOP)/targetlib.mk


