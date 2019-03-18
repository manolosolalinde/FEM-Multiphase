# $Id: Makefile 3101 2008-10-16 19:34:03Z roystgnr $


# The location of the mesh library
# LIBMESH_DIR ?= /opt/libmesh/libmesh_4180
LIBMESH_DIR ?= /mnt/0F22134B0F22134B/GITHUB/petsc/arch-linux2-c-debug/externalpackages/git.libmesh
BUBBLEPROJECT_DIR ?= ./
#LIBMESH_RUN = mpirun -np 4
#mpirun -np 2 ./ex19-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
LIBMESH_OPTIONS = -ksp_converged_reason #-ksp_type preonly -pc_type redistribute -redistribute_ksp_type gmres #-redistribute_pc_type asm -log_summary
#-snes_converged_reason -snes_mf_operator

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(LIBMESH_DIR)/Make.common
#include /opt/netcdf/netcdf-4.1.1/cxx


###############################################################################
# File management.  This is where the source, header, and object files are
# defined
out1 = testreinit.C test.C testadvection.C
out2 = main.C testreinit.C testadvection.C
out3 = test.C main.C testadvection.C
out4 = test.C main.C testreinit.C
#
# source files
srcfiles_main	:=  $(filter-out $(out1),$(wildcard *.C))
srcfiles_2      :=  $(filter-out $(out2),$(wildcard *.C))
srcfiles_3      :=  $(filter-out $(out3),$(wildcard *.C))
srcfiles_4      :=  $(filter-out $(out4),$(wildcard *.C))

#
# object files
objects_main	:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles_main))
objects_2	:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles_2))
objects_3	:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles_3))
objects_4	:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles_4))
###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
backup_target := ./backups/backuper
sim1_target   := ./sim1-$(METHOD)
sim2_target   := ./sim2-$(METHOD)
sim3_target   := ./sim3-$(METHOD)
sim4_target   := ./sim4-$(METHOD)
testreinit_target := ./testreinit-$(METHOD)
testadvection_target := ./testadvection-$(METHOD)


# Production rules:  how to make the target - depends on library configuration
$(sim1_target): $(objects_main)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects_main) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)

$(sim2_target): $(objects_main)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects_main) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)

$(sim3_target): $(objects_main)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects_main) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)

$(sim4_target): $(objects_main)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects_main) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)

$(testreinit_target): $(objects_3)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects_3) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)

$(testadvection_target): $(objects_4)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects_4) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)

# Useful rules.
clean:
	@rm -f *.o *.g.o *.pg.o *~ .depend

cleanexe:
	@rm -f *-opt *-dbg

distclean:
	@$(MAKE) clean
	@$(MAKE) cleanexe

backup: 
	@ $(backup_target)

# Opciones PETSC -ksp_type preonly -pc_type lu -pc_factor_shift_nonzero

sim1: $(sim1_target)
	@echo "***************************************************************"
	@echo "* Running Multiphase Test " $(LIBMESH_RUN) $(sim1_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESH_RUN) $(sim1_target) $(LIBMESH_OPTIONS) | tee outfile_sim1.txt
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Multiphase Test " $(LIBMESH_RUN) $(sim1_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"

sim2: $(sim2_target)
	@echo "***************************************************************"
	@echo "* Running Multiphase Test " $(LIBMESH_RUN) $(sim2_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESH_RUN) $(sim2_target) $(LIBMESH_OPTIONS) | tee outfile_sim2.txt
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Multiphase Test " $(LIBMESH_RUN) $(sim2_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"

sim3: $(sim3_target)
	@echo "***************************************************************"
	@echo "* Running Multiphase Test " $(LIBMESH_RUN) $(sim3_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESH_RUN) $(sim3_target) $(LIBMESH_OPTIONS) | tee outfile_sim3.txt
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Multiphase Test " $(LIBMESH_RUN) $(sim3_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"

sim4: $(sim4_target)
	@echo "***************************************************************"
	@echo "* Running Multiphase Test " $(LIBMESH_RUN) $(sim4_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESH_RUN) $(sim4_target) $(LIBMESH_OPTIONS) | tee outfile_sim4.txt
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Multiphase Test " $(LIBMESH_RUN) $(sim4_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"

testreinit: $(testreinit_target)
	@echo "***************************************************************"
	@echo "* Running Reinitialization Test " $(LIBMESH_RUN) $(testreinit_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESH_RUN) $(testreinit_target) $(LIBMESH_OPTIONS) --enable-reference-counting | tee outfile_reinit.txt
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Reinitialization Test" $(LIBMESH_RUN) $(testreinit_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"

testadvection: $(testadvection_target)
	@echo "***************************************************************"
	@echo "* Running Advection Test" $(LIBMESH_RUN) $(testadvection_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESH_RUN) $(testadvection_target) $(LIBMESH_OPTIONS) -ksp_type preonly -pc_type lu -pc_factor_shift_nonzero --enable-reference-counting
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Advection Test " $(LIBMESH_RUN) $(testadvection_target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"


# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@$(perl) -pi -e 's#    $(LIBMESH_DIR)#    \$$\(LIBMESH_DIR\)#' .depend
	@echo "Updated .depend"

###############################################################################
