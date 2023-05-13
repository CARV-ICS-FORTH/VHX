#
# Copyright (c) 2021-2022 Computer Architecture and VLSI Systems (CARV)
#                         Laboratory, ICS Forth. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#


sources = \
	coll_vhx.h \
	coll_vhx.c \
	coll_vhx_component.c \
	coll_vhx_module.c \
	coll_vhx_barrier.c\
	coll_vhx_bcast.c \
	coll_vhx_allreduce.c \
	coll_vhx_reduce.c
	
# Make the output library in this directory, and name it either
# mca_<type>_<name>.la (for DSO builds) or libmca_<type>_<name>.la
# (for static builds).

component_noinst =
component_install =
if MCA_BUILD_ompi_coll_vhx_DSO
component_install += mca_coll_vhx.la
else
component_noinst += libmca_coll_vhx.la
endif

mcacomponentdir = $(ompilibdir)
mcacomponent_LTLIBRARIES = $(component_install)
mca_coll_vhx_la_SOURCES = $(sources)
mca_coll_vhx_la_LDFLAGS = -module -avoid-version
mca_coll_vhx_la_LIBADD = $(top_builddir)/ompi/lib@OMPI_LIBMPI_NAME@.la

noinst_LTLIBRARIES = $(component_noinst)
libmca_coll_vhx_la_SOURCES = $(sources)
libmca_coll_vhx_la_LDFLAGS = -module -avoid-version
