AMREX_HOME ?= ../../amrex
KOKKOS_HOME ?= ../../amrex/Tools/kokkoslib/install

DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_CUDA  = TRUE

USE_MPI   = FALSE
USE_OMP   = FALSE

TINY_PROFILE = TRUE

INCLUDE_LOCATIONS += $(KOKKOS_HOME)/include
LIBRARY_LOCATIONS += $(KOKKOS_HOME)/lib 
LIBRARIES += -lkokkos

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
