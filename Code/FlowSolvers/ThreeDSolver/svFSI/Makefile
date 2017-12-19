# Hey emacs, this is a -*- makefile -*-

# Copyright (c) 2014-2015 The Regents of the University of California.
# All Rights Reserved. 
#
# Portions of the code Copyright (c) 2009-2011 Open Source Medical
# Software Corporation, University of California, San Diego. 
#
# Portions Copyright (c) 2000-2007, Stanford University, 
#    Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
#    Charles A. Taylor.
#
# See SimVascular Acknowledgements file for additional
# contributors to the source code.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions 
# are met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer. 
# Redistributions in binary form must reproduce the above copyright 
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution. 
# Neither the name of the Stanford University or Rensselaer Polytechnic
# Institute nor the names of its contributors may be used to endorse or
# promote products derived from this software without specific prior 
# written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

TOP=../../../../BuildWithMake

include $(TOP)/include.mk

# --------------
# Compiler flags
# --------------

CXXFLAGS = $(GLOBAL_CXXFLAGS) $(MPI_INCDIR) $(PARMETIS_INCDIR) $(METIS_INCDIR) $(LESLIB_INCDIR) $(LESLIB_DEFS) $(SVLS_INCDIR) $(SVLS_DEFS) $(SOLVERIO_INCDIR) $(THREEDSOLVER_INCDIR) $(ZLIB_INCDIR) $(VTK_INCDIRS) -I../../Include/Make
CCFLAGS  = $(GLOBAL_CCFLAGS) $(MPI_INCDIR) $(PARMETIS_INCDIR) $(METIS_INCDIR) $(LESLIB_INCDIR) $(LESLIB_DEFS) $(SVLS_INCDIR) $(SVLS_DEFS) $(SOLVERIO_INCDIR) $(THREEDSOLVER_INCDIR) $(ZLIB_INCDIR) $(VTK_INCDIRS) -I../../Include/Make
FFLAGS   = $(GLOBAL_FFLAGS) $(MPI_INCDIR) $(PARMETIS_INCDIR) $(METIS_INCDIR) $(SVLS_INCDIR) $(SVLS_DEFS) $(THREEDSOLVER_INCDIR) -I../../Include/Make

# c headers here
HDRS	= 

#CXXSRCS = remeshTet.cxx remeshSim.cxx

# for shape functions
CSRCS    = vtkZpipe.c SPLIT.c

# platform specific functions
#ifneq ($(CLUSTER), x64_linux)
#  CSRCS   += missing_win32.c 
#endif

# must build modules first!

FSRCS   = UTIL.f \
         COMU.f \
         CHNL.f \
         LIST.f \
         MOD.f \
         MATFUN.f \
         ALLFUN.f \
         BBO.f \
         BAFINI.f \
         CONSTRUCT.f \
         DISTRIBUTE.f \
         FFT.f \
         FLUID.f \
         HEAT.f \
         INITIALIZE.f \
         LELAS.f \
         LHSA.f \
         LOADMSH.f \
         LOADNRB.f \
         MAIN.f \
         MATMODELS.f \
         NN.f \
         NURBS.f \
         OUTPUT.f \
         PIC.f \
         POST.f \
         READFILES.f \
         READMSH.f \
         REMESH.f \
         SETBC.f \
         STRUCT.f \
         TXT.f \
         VTKXML.f

F90SRCS = vtkComMod.f90 \
          vtkTypeParams.f90 \
          vtkXMLParser.f90  \
          vtkLegacyParser.f90

SRCS	= $(FSRCS) $(CSRCS) $(CXXSRCS)

DEPS	= $(CSRCS:.c=.d) $(CXXSRCS:.cxx=.d)

OBJS    = $(addprefix $(BUILD_MPI_DIR)/,$(CXXSRCS:.cxx=.$(OBJECTEXT))) \
          $(addprefix $(BUILD_MPI_DIR)/,$(CSRCS:.c=.$(OBJECTEXT))) \
          $(addprefix $(BUILD_MPI_DIR)/,$(F90SRCS:.f90=.$(OBJECTEXT))) \
          $(addprefix $(BUILD_MPI_DIR)/,$(FSRCS:.f=.$(OBJECTEXT)))

TARGET_EXE_FULL = $(TOP)/Bin/svfsi-$(CXX_COMPILER_VERSION)-$(FORTRAN_COMPILER_VERSION)-$(MPI_NAME).exe
TARGET_EXE = $(TOP)/Bin/svfsi-$(MPI_NAME).exe

static: directories $(TARGET_EXE_FULL) 

directories:
	-mkdir -p $(BUILD_MPI_DIR)

$(TARGET_EXE_FULL): $(OBJS)
	for fn in $(TARGET_EXE); do /bin/rm -f $$fn; done
	for fn in $(basename $(TARGET_EXE)).pdb; do /bin/rm -f $$fn; done
	for fn in $(TARGET_EXE_FULL); do /bin/rm -f $$fn; done
	for fn in $(basename $(TARGET_EXE_FULL)).pdb; do /bin/rm -f $$fn; done
	$(LINK_EXE)$(TARGET_EXE_FULL) $(GLOBAL_LFLAGS)  $(SOLVER_STCKSZ) \
             $(OBJS) \
             $(SVLIBFLAG)_simvascular_solverio$(LIBLINKEXT) \
	     $(ZLIB_LIBS) $(VTK_LIBS) $(PARMETIS_LIBS) $(METIS_LIBS) $(LESLIB_LIBS) \
             $(SVLS_LIBS) $(MPI_LIBS) $(F90_LIBS) $(CXX_LIBS)
	if [ -e $(notdir $(TARGET_EXE_FULL)).pdb ];then mv $(basename $(TARGET_EXE_FULL)).pdb $(TOP)/Bin;fi
	ln -s $(notdir $(TARGET_EXE_FULL)) $(TARGET_EXE)

ifndef NO_DEPEND
-include $(DEPS)
endif

clean:
	for fn in $(BUILD_DIR); do /bin/rm -f -r $$fn;done
	for fn in *~; do /bin/rm -f $$fn;done
	for fn in *.mod; do /bin/rm -f $$fn;done
	for fn in $(TOP)/Bin/$(TARGET_EXE)*; do /bin/rm -f $$fn; done

veryclean: clean
	if [ -e obj ];then /bin/rm -f -r obj;fi

