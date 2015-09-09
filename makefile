Here = $(PWD)
ConfigFile = $(Here)/ttbarjets.cfg
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
DipoleDir = $(Here)/dipoles
PDFDir = $(Here)/PDFS
VegasDir = $(Here)/Vegas
JHUGenMELADir = $(Here)/JHUGenMELA
OptReport = $(Here)/OptRep.txt
PSDir = $(Here)/PhaseSpace
QCDLoop = $(Here)/QCDLoop-1.9
Flags = 

# ifort optimization, Yes/No
Opt = Yes

# MPI features, Yes/No. run with: mpiexec -n 4 ./TOPAZ_MPI ...
useMPI = No

# link pdfs via LHA library ('Yes' or 'No')
UseLHAPDF=No
LHAPDFDir=/afs/cern.ch/user/m/maschulz/lib/LHAPDF-6.1.5/lib/
# LHAPDFDir=directory which contains libLHAPDF.a, libLHAPDF.la, libLHAPDF.so
# remember to export 
#          LD_LIBRARY_PATH=/.../LHAPDF-x.y.z/lib/:${LD_LIBRARY_PATH}
#          LHAPDF_DATA_PATH=/.../LHAPDF-x.y.z/share/LHAPDF/:${LHAPDF_DATA_PATH}


# interface to JHUGenMELA (requires UseLHAPDF=Yes)
useJHUGenMELA = No


ifeq ($(useMPI),Yes)
    Exec = ./TOPAZ_MPI
    F95compiler = mpif90 -f90=ifort -lpthread -D_UseMPIVegas=1
    ccomp = mpicc -lpthread  -lm 
else
    Exec = ./TOPAZ
    F95compiler = ifort -D_UseMPIVegas=0
    ccomp = gcc -O2
endif





ifeq ($(UseLHAPDF),Yes)
    LHAPDFflags = -L$(LHAPDFDir) -lLHAPDF -D_UseLHAPDF=1
else
    LHAPDFflags = -D_UseLHAPDF=0
endif 
 
 
ifeq ($(useJHUGenMELA),Yes)
    Flags += -D_UseJHUGenMELA=1
else
    Flags += -D_UseJHUGenMELA=0
endif
 

ifeq ($(Opt),Yes)
   IfortOpts   = -O2 -fpp -opt-report -opt-report-file$(OptReport) -I$(Here)/colors -I$(VegasDir) -module $(ModuleDir) $(LHAPDFflags) $(Flags)
else
   IfortOpts   = -O0 -fpp -implicitnone -check bounds -check pointer -warn interfaces -ftrapuv  -debug extended -g -traceback -fpe0 -check uninit -I$(Here)/colors -I$(VegasDir) -module $(ModuleDir) $(LHAPDFflags) $(Flags)
endif
fcomp = $(F95compiler) $(IfortOpts) @$(ConfigFile)






# makeDep = $(ConfigFile) makefile


# fastjet stuff
#FASTJET_CONFIG=/home/schulze/usr/local/bin/fastjet-config
# CXXFLAGS += $(shell $(FASTJET_CONFIG) --cxxflags)
# FJLIBS += $(shell $(FASTJET_CONFIG) --libs --plugins )


PDFObj = $(ObjectDir)/mstwpdf.o \
         $(ObjectDir)/cteq2mrst.o \
         $(ObjectDir)/mrst2001lo.o \
         $(ObjectDir)/Cteq66Pdf.o \
         $(ObjectDir)/CT10Pdf.o \
         $(ObjectDir)/NNPDFDriver.o
         
RockyObj = $(ObjectDir)/genps.o \
           $(ObjectDir)/boost.o
           
YetiObj  = $(ObjectDir)/yeti.o

ifeq ($(useMPI),Yes)
   VegasObj = $(VegasDir)/pvegas_mpi.o
else
   VegasObj = $(VegasDir)/vegas.o 
endif


# JHUGenMELA modules 
ifeq ($(useJHUGenMELA),Yes)
   JHUGenMELADep = $(JHUGenMELADir)/mod_TTBH_MatEl.F90
   JHUGenMELAObj = $(ObjectDir)/mod_TTBH_MatEl.o           
else
   JHUGenMELADep = 
   JHUGenMELAObj = 
endif

IntegralObj = $(QCDLoop)/ql/libqcdloop.a\
              $(QCDLoop)/ff/libff.a


OPPDep = mod_NVBasis.f90 \
         mod_Residues.f90 \
         mod_UCuts.f90 \
         mod_Residues_new.f90 \
         mod_UCuts_new.f90

OPPObj = $(ObjectDir)/mod_NVBasis.o \
         $(ObjectDir)/mod_Residues.o \
         $(ObjectDir)/mod_UCuts.o \
         $(ObjectDir)/mod_Residues_new.o \
         $(ObjectDir)/mod_UCuts_new.o \
         $(ObjectDir)/mod_NVBasis128.o \
         $(ObjectDir)/mod_Residues128.o \
         $(ObjectDir)/mod_UCuts128.o \
         $(ObjectDir)/mod_Residues128_new.o \
         $(ObjectDir)/mod_UCuts128_new.o 

DipoleDepTTB = $(Here)/dipoles/mod_Dipoles_GGTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_GGTTBG_noDK.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQ.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQ_noDK.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBG_noDK.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBG_noDK.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQ.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQ_noDK.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBG_noDK.f90

DipoleObjTTB = $(ObjectDir)/mod_Dipoles_GGTTBG.o \
            $(ObjectDir)/mod_Dipoles_GGTTBG_noDK.o \
            $(ObjectDir)/mod_Dipoles_QGTTBQ.o \
            $(ObjectDir)/mod_Dipoles_QGTTBQ_noDK.o \
            $(ObjectDir)/mod_Dipoles_QQBTTBG.o \
            $(ObjectDir)/mod_Dipoles_QQBTTBG_noDK.o \
            $(ObjectDir)/mod_IntDipoles_GGTTBG.o \
            $(ObjectDir)/mod_IntDipoles_GGTTBG_noDK.o \
            $(ObjectDir)/mod_IntDipoles_QGTTBQ.o \
            $(ObjectDir)/mod_IntDipoles_QGTTBQ_noDK.o \
            $(ObjectDir)/mod_IntDipoles_QQBTTBG.o \
            $(ObjectDir)/mod_IntDipoles_QQBTTBG_noDK.o




DipoleDepTTBJ = $(Here)/dipoles/mod_Dipoles_GGTTBGG2.f90 \
            $(Here)/dipoles/mod_Dipoles_GGTTBQQB2.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBGG2.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQG2.f90 \
            $(Here)/dipoles/mod_Dipoles_QBGTTBQBG2.f90 \
            $(Here)/dipoles/mod_Dipoles_qbb_ttqqqq_dip2.f90 \
            $(Here)/dipoles/mod_Dipoles_qqb_ttqqqq_dip2.f90 \
            $(Here)/dipoles/mod_Dipoles_qqq_ttqqqq_dip2.f90 \
            $(Here)/dipoles/mod_Dipoles_DKJ_TTB.f90 \
            $(Here)/dipoles/mod_Dipoles_DKJ_GGTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_DKJ_QQBTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_DKJ_QGTTBQ.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBGG2.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBQQB2.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBGG2.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQG2.f90 \
            $(Here)/dipoles/mod_IntDipoles_QBGTTBQBG2.f90 \
            $(Here)/dipoles/mod_IntDipoles_SixQuark2.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKJ_GGTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKJ_QQBTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKJ_QGTTBQ.f90


DipoleObjTTBJ = $(ObjectDir)/mod_Dipoles_GGTTBGG2.o \
            $(ObjectDir)/mod_Dipoles_GGTTBQQB2.o \
            $(ObjectDir)/mod_Dipoles_QQBTTBGG2.o \
            $(ObjectDir)/mod_Dipoles_QGTTBQG2.o \
            $(ObjectDir)/mod_Dipoles_QBGTTBQBG2.o \
            $(ObjectDir)/mod_Dipoles_qbb_ttqqqq_dip2.o \
            $(ObjectDir)/mod_Dipoles_qqb_ttqqqq_dip2.o \
            $(ObjectDir)/mod_Dipoles_qqq_ttqqqq_dip2.o \
            $(ObjectDir)/mod_Dipoles_DKJ_TTB.o \
            $(ObjectDir)/mod_Dipoles_DKJ_GGTTBG.o \
            $(ObjectDir)/mod_Dipoles_DKJ_QQBTTBG.o \
            $(ObjectDir)/mod_Dipoles_DKJ_QGTTBQ.o \
            $(ObjectDir)/mod_IntDipoles_GGTTBGG2.o \
            $(ObjectDir)/mod_IntDipoles_GGTTBQQB2.o \
            $(ObjectDir)/mod_IntDipoles_QQBTTBGG2.o \
            $(ObjectDir)/mod_IntDipoles_QGTTBQG2.o \
            $(ObjectDir)/mod_IntDipoles_QBGTTBQBG2.o \
            $(ObjectDir)/mod_IntDipoles_SixQuark2.o \
            $(ObjectDir)/mod_IntDipoles_DKJ_GGTTBG.o \
            $(ObjectDir)/mod_IntDipoles_DKJ_QQBTTBG.o \
            $(ObjectDir)/mod_IntDipoles_DKJ_QGTTBQ.o





DipoleDepTTBP = $(Here)/dipoles/mod_Dipoles_GGTTBGP.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBGP.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQP.f90 \
            $(Here)/dipoles/mod_Dipoles_QBGTTBQBP.f90 \
            $(Here)/dipoles/mod_Dipoles_DKP_GGTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_DKP_QQBTTBG.f90 \
            $(Here)/dipoles/mod_Dipoles_DKP_QGTTBQ.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBGP.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBGP.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQP.f90 \
            $(Here)/dipoles/mod_IntDipoles_QBGTTBQBP.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKP_GGTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKP_QQBTTBG.f90 \
            $(Here)/dipoles/mod_IntDipoles_DKP_QGTTBQ.f90

DipoleObjTTBP = $(ObjectDir)/mod_Dipoles_GGTTBGP.o \
                $(ObjectDir)/mod_Dipoles_QQBTTBGP.o \
                $(ObjectDir)/mod_Dipoles_QGTTBQP.o \
                $(ObjectDir)/mod_Dipoles_QBGTTBQBP.o \
                $(ObjectDir)/mod_Dipoles_DKP_GGTTBG.o \
                $(ObjectDir)/mod_Dipoles_DKP_QQBTTBG.o \
                $(ObjectDir)/mod_Dipoles_DKP_QGTTBQ.o \
                $(ObjectDir)/mod_IntDipoles_GGTTBGP.o \
                $(ObjectDir)/mod_IntDipoles_QQBTTBGP.o \
                $(ObjectDir)/mod_IntDipoles_QGTTBQP.o \
                $(ObjectDir)/mod_IntDipoles_QBGTTBQBP.o \
                $(ObjectDir)/mod_IntDipoles_DKP_GGTTBG.o \
                $(ObjectDir)/mod_IntDipoles_DKP_QQBTTBG.o \
                $(ObjectDir)/mod_IntDipoles_DKP_QGTTBQ.o



DipoleDepTTBZ = $(Here)/dipoles/mod_Dipoles_GGTTBGZ.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBGZ.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQZ.f90 \
            $(Here)/dipoles/mod_Dipoles_QBGTTBQBZ.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBGZ.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBGZ.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQZ.f90 \
            $(Here)/dipoles/mod_IntDipoles_QBGTTBQBZ.f90

DipoleObjTTBZ = $(ObjectDir)/mod_Dipoles_GGTTBGZ.o \
                $(ObjectDir)/mod_Dipoles_QQBTTBGZ.o \
                $(ObjectDir)/mod_Dipoles_QGTTBQZ.o \
                $(ObjectDir)/mod_Dipoles_QBGTTBQBZ.o \
                $(ObjectDir)/mod_IntDipoles_GGTTBGZ.o \
                $(ObjectDir)/mod_IntDipoles_QQBTTBGZ.o \
                $(ObjectDir)/mod_IntDipoles_QGTTBQZ.o \
                $(ObjectDir)/mod_IntDipoles_QBGTTBQBZ.o


DipoleDepTTBH = $(Here)/dipoles/mod_Dipoles_GGTTBGH.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBGH.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQH.f90 \
            $(Here)/dipoles/mod_Dipoles_QBGTTBQBH.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBGH.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBGH.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQH.f90 \
            $(Here)/dipoles/mod_IntDipoles_QBGTTBQBH.f90

DipoleObjTTBH = $(ObjectDir)/mod_Dipoles_GGTTBGH.o \
                $(ObjectDir)/mod_Dipoles_QQBTTBGH.o \
                $(ObjectDir)/mod_Dipoles_QGTTBQH.o \
                $(ObjectDir)/mod_Dipoles_QBGTTBQBH.o \
                $(ObjectDir)/mod_IntDipoles_GGTTBGH.o \
                $(ObjectDir)/mod_IntDipoles_QQBTTBGH.o \
                $(ObjectDir)/mod_IntDipoles_QGTTBQH.o \
                $(ObjectDir)/mod_IntDipoles_QBGTTBQBH.o

DipoleDepTTBP_anom = $(Here)/dipoles/mod_Dipoles_GGTTBGP_anom.f90 \
            $(Here)/dipoles/mod_Dipoles_QQBTTBGP_anom.f90 \
            $(Here)/dipoles/mod_Dipoles_QGTTBQP_anom.f90 \
            $(Here)/dipoles/mod_Dipoles_QBGTTBQBP_anom.f90 \
            $(Here)/dipoles/mod_IntDipoles_GGTTBGP_anom.f90 \
            $(Here)/dipoles/mod_IntDipoles_QQBTTBGP_anom.f90 \
            $(Here)/dipoles/mod_IntDipoles_QGTTBQP_anom.f90 \
            $(Here)/dipoles/mod_IntDipoles_QBGTTBQBP_anom.f90

DipoleObjTTBP_anom = $(ObjectDir)/mod_Dipoles_GGTTBGP_anom.o \
                $(ObjectDir)/mod_Dipoles_QQBTTBGP_anom.o \
                $(ObjectDir)/mod_Dipoles_QGTTBQP_anom.o \
                $(ObjectDir)/mod_Dipoles_QBGTTBQBP_anom.o \
                $(ObjectDir)/mod_IntDipoles_GGTTBGP_anom.o \
                $(ObjectDir)/mod_IntDipoles_QQBTTBGP_anom.o \
                $(ObjectDir)/mod_IntDipoles_QGTTBQP_anom.o \
                $(ObjectDir)/mod_IntDipoles_QBGTTBQBP_anom.o

                
DipoleDepSTSTB = $(Here)/dipoles/mod_Dipoles_GGSTSTBG.f90 \
		 $(Here)/dipoles/mod_IntDipoles_GGSTSTBG.f90 \
		 $(Here)/dipoles/mod_Dipoles_QQBSTSTBG.f90  \
		 $(Here)/dipoles/mod_Dipoles_QGSTSTBQ.f90 \
                 $(Here)/dipoles/mod_IntDipoles_QQBSTSTBG.f90 \
                 $(Here)/dipoles/mod_IntDipoles_QGSTSTBQ.f90

DipoleObjSTSTB = $(ObjectDir)/mod_Dipoles_GGSTSTBG.o \
		 $(ObjectDir)/mod_IntDipoles_GGSTSTBG.o \
		 $(ObjectDir)/mod_Dipoles_QQBSTSTBG.o  \
		 $(ObjectDir)/mod_Dipoles_QGSTSTBQ.o \
                 $(ObjectDir)/mod_IntDipoles_QQBSTSTBG.o \
                 $(ObjectDir)/mod_IntDipoles_QGSTSTBQ.o




DipoleDepHTHTB = $(Here)/dipoles/mod_Dipoles_GGHTHTBG.f90 \
		 $(Here)/dipoles/mod_Dipoles_QQBHTHTBG.f90  \
		 $(Here)/dipoles/mod_Dipoles_QGHTHTBQ.f90  \
		 $(Here)/dipoles/mod_IntDipoles_GGHTHTBG.f90 \
                 $(Here)/dipoles/mod_IntDipoles_QQBHTHTBG.f90 \
                 $(Here)/dipoles/mod_IntDipoles_QGHTHTBQ.f90

DipoleObjHTHTB = $(ObjectDir)/mod_Dipoles_GGHTHTBG.o \
		 $(ObjectDir)/mod_Dipoles_QQBHTHTBG.o  \
		 $(ObjectDir)/mod_Dipoles_QGHTHTBQ.o \
		 $(ObjectDir)/mod_IntDipoles_GGHTHTBG.o \
                 $(ObjectDir)/mod_IntDipoles_QQBHTHTBG.o \
                 $(ObjectDir)/mod_IntDipoles_QGHTHTBQ.o



DipoleDepZprime = $(Here)/dipoles/mod_Dipoles_ZprimeTTB.f90 \
		  $(Here)/dipoles/mod_IntDipoles_ZprimeTTB.f90

DipoleObjZprime = $(ObjectDir)/mod_Dipoles_ZprimeTTB.o \
		 $(ObjectDir)/mod_IntDipoles_ZprimeTTB.o 




DipoleDepeeTTB = $(Here)/dipoles/mod_Dipoles_eeTTB.f90 \
		  $(Here)/dipoles/mod_IntDipoles_eeTTB.f90

DipoleObjeeTTB = $(ObjectDir)/mod_Dipoles_eeTTB.o \
		 $(ObjectDir)/mod_IntDipoles_eeTTB.o 

		 
		 
MadGraphObj = $(Here)/MadGraph/gg_ttbg.o \
				  $(Here)/MadGraph/switchmom.o \
				  $(HOME)/lib/HELAS-3.0/coupsm.o \
				  $(HOME)/lib/HELAS-3.0/oxxxxx.o \
				  $(HOME)/lib/HELAS-3.0/vxxxxx.o \
				  $(HOME)/lib/HELAS-3.0/ixxxxx.o \
				  $(HOME)/lib/HELAS-3.0/fvixxx.o \
				  $(HOME)/lib/HELAS-3.0/sxxxxx.o \
				  $(HOME)/lib/HELAS-3.0/fsoxxx.o \
				  $(HOME)/lib/HELAS-3.0/iosxxx.o \
				  $(HOME)/lib/HELAS-3.0/fsixxx.o \
				  $(HOME)/lib/HELAS-3.0/jvvxxx.o \
				  $(HOME)/lib/HELAS-3.0/jgggxx.o \
				  $(HOME)/lib/HELAS-3.0/jioxxx.o \
				  $(HOME)/lib/HELAS-3.0/iovxxx.o \
				  $(HOME)/lib/HELAS-3.0/fvoxxx.o \
				  $(HOME)/lib/HELAS-3.0/vvvxxx.o \
				  $(HOME)/lib/HELAS-3.0/libdhelas3.ifc90.a \
				  $(Here)/MadGraph/gg_ttbh.o\
				  $(Here)/MadGraph/uub_ttbh.o\
				  $(Here)/MadGraph/gg_ttbhg.o\
				  $(Here)/MadGraph/uub_ttbhg.o\
				  $(Here)/MadGraph/ubg_ttbhub.o \
				  $(Here)/MadGraph/ug_ttbhu.o
# 				  $(Here)/MadGraph/emep_ttb.o \
# 				  $(Here)/MadGraph/emep_ttbg.o \
# 				  $(Here)/MadGraph/gg_ttbz.o \
# 				  $(Here)/MadGraph/gg_ttbzg.o\
# 				  $(Here)/MadGraph/uub_ttbgz.o\
# 				  $(Here)/MadGraph/ddb_ttbgz.o \
# 				  $(Here)/MadGraph/dg_ttbdz.o \
# 				  $(Here)/MadGraph/gub_ttbubz.o \
# 				  $(Here)/MadGraph/gdb_ttbdbz.o \
# 				  $(Here)/MadGraph/ug_ttbuz.o \
# 				  $(Here)/MadGraph/uub_ttbz.o \
# 				  $(Here)/MadGraph/ddb_ttbz.o \
# 				  $(Here)/MadGraph/ubg_ttbubz.o
#				  $(Here)/MadGraph/gg_ttb.o \
#				  $(Here)/MadGraph/gg_ttba.o \
#				  $(Here)/MadGraph/gg_tbtga.o \
#				  $(Here)/MadGraph/ug_ttbua.o \
#				  $(Here)/MadGraph/ddb_ttba.o \
#				  $(Here)/MadGraph/uub_ttba.o \
#				  $(Here)/MadGraph/dbg_ttbdba.o \
#				  $(Here)/MadGraph/ddb_ttbga.o \
#				  $(Here)/MadGraph/wm_ubduub.o  \
#				  $(Here)/MadGraph/wm_ubdccb.o  \
#				  $(Here)/MadGraph/gg_ttbgg.o \
#				  $(Here)/MadGraph/tb_bbemveb.o \
#				  $(Here)/MadGraph/tb_tbepem.o \
#				  $(Here)/MadGraph/gg_bepvebbemve.o \
#				  $(Here)/MadGraph/gg_uub.o \
#				  $(Here)/MadGraph/gg_uubg.o \
#				  $(Here)/MadGraph/uub_ttb.o \
#				  $(Here)/MadGraph/uub_ttbg.o \
#				  $(Here)/MadGraph/uub_ttbgg.o \
#				  $(Here)/MadGraph/uub_ttbddb.o \
#				  $(Here)/MadGraph/uub_ttbuub.o \
#				  $(Here)/MadGraph/ubdb_ttbubdb.o \
#				  $(Here)/MadGraph/udb_ttbudb.o \
#				  $(Here)/MadGraph/ud_ttbud.o \
#				  $(Here)/MadGraph/uu_ttbuu.o \
#				  $(Here)/MadGraph/ubub_ttbubub.o \
#				  $(Here)/MadGraph/gg_ttbuub.o \
#				  $(Here)/MadGraph/uub_ddbg.o \
#				  $(Here)/MadGraph/ug_ttbu.o \
#				  $(Here)/MadGraph/ubg_ttbub.o \
#				  $(Here)/MadGraph/gub_ttbub.o \
#				  $(Here)/MadGraph/ug_ttbug.o \
#				  $(Here)/MadGraph/gu_ttbu.o \
#				  $(Here)/MadGraph/ubg_ttbubg.o \
#				  $(Here)/MadGraph/t_bepve.o \
#				  $(Here)/MadGraph/tb_emvebbb.o \
#				  $(Here)/MadGraph/tb_bbemveb.o \
#				  $(Here)/MadGraph/tb_bbemvebgg.o \
#				  $(Here)/MadGraph/t_bepvegg.o \
#				  $(Here)/MadGraph/t_epvebbbb.o \
#				  $(Here)/MadGraph/t_epvebuub.o \
#				  $(Here)/MadGraph/tb_emvebbbuub.o \
#				  $(Here)/MadGraph/tb_emvebbbbbb.o \
#				  $(Here)/MadGraph/t_bdbuccb.o \
#				  $(Here)/MadGraph/t_bdbuddb.o \
#				  $(Here)/MadGraph/t_bdbuuub.o \
#				  $(Here)/MadGraph/tb_bbubdccb.o \
#				  $(Here)/MadGraph/tb_bbubdddb.o \
#				  $(Here)/MadGraph/tb_bbubduub.o \
#				  $(Here)/MadGraph/tb_bbubdgg.o \
#				  $(Here)/MadGraph/t_bdbugg.o \
#				  $(Here)/MadGraph/tb_bbwm.o  \
#				  $(Here)/MadGraph/t_bwp.o \
#				  $(Here)/MadGraph/t_epveb.o \
#				  $(Here)/MadGraph/t_bwpuub.o \
#				  $(Here)/MadGraph/t_bwpbbb.o \
#				  $(Here)/MadGraph/wm_emveb.o \
#				  $(Here)/MadGraph/wp_epve.o \


# ------------------------------------------------------------


# the order of these object files corresponds to their mutual dependency
allObjects =   				$(ObjectDir)/mod_Misc.o \
					$(ObjectDir)/mod_Parameters.o \
					$(ObjectDir)/mod_Process.o \
					$(ObjectDir)/mod_Permutations.o \
					$(ObjectDir)/mod_IntegerPartition.o \
					$(ObjectDir)/mod_MyRecurrence.o \
					$(ObjectDir)/mod_MyWeylRecurrence.o \
					$(ObjectDir)/mod_Amplitudes.o \
					$(ObjectDir)/mod_SingleTopHAmps.o \
					$(ObjectDir)/mod_Amplitudes_Zprime.o \
					$(ObjectDir)/mod_Amplitudes_eeTTB.o \
					$(OPPObj) \
					$(ObjectDir)/mod_DKIntDipoles.o \
					$(ObjectDir)/mod_HadrWDecay.o \
					$(ObjectDir)/mod_JPsiFrag.o \
					$(Here)/includes_DKP/DKP1L.o \
					$(Here)/includes_DKJ/TopCoeffsNLO.o \
					$(Here)/includes_DKJ/ATopCoeffsNLO.o \
					$(ObjectDir)/mod_TTBP_NLODK.o \
					$(ObjectDir)/mod_TTBJ_NLODK.o \
					$(ObjectDir)/mod_TTBJ_NLODKW.o \
					$(ObjectDir)/mod_WDecay.o \
					$(ObjectDir)/mod_ZDecay.o \
					$(ObjectDir)/mod_HDecay.o \
					$(ObjectDir)/mod_TopDecay.o \
					$(ObjectDir)/mod_ExoticDecay.o \
					$(ObjectDir)/mod_IntDipoles.o \
					$(ObjectDir)/mod_SpinCorrel.o \
					$(ObjectDir)/mod_Kinematics.o \
					$(ObjectDir)/mod_Weighting.o \
					$(ObjectDir)/mod_Integrals.o \
					$(DipoleObjTTB) \
					$(DipoleObjTTBJ) \
					$(DipoleObjTTBP) \
					$(DipoleObjTTBZ) \
					$(DipoleObjTTBH) \
					$(DipoleObjTTBP_anom) \
					$(DipoleObjSTSTB) \
					$(DipoleObjHTHTB) \
					$(DipoleObjZprime) \
					$(DipoleObjeeTTB) \
					$(ObjectDir)/mod_SixFermionProcs2.o \
					$(ObjectDir)/mod_CrossSection_TTB.o \
					$(ObjectDir)/mod_CrossSection_TTBJ.o \
					$(ObjectDir)/mod_CrossSection_TTBP.o \
					$(ObjectDir)/mod_CrossSection_TTBP_anomcoupl.o \
					$(ObjectDir)/mod_CrossSection_TTBETmiss.o \
					$(ObjectDir)/mod_CrossSection_ZprimeTTB.o \
					$(ObjectDir)/mod_CrossSection_eeTTB.o \
					$(ObjectDir)/mod_CrossSection_TTBZ.o \
					$(ObjectDir)/mod_CrossSection_TTBH.o \
					$(ObjectDir)/mod_CrossSection_TH.o \
					$(ObjectDir)/main.o


#--------------------------------------------------------------------------------------------------
# note that DKP1L.o,TopCoeffsNLO.o,ATopCoeffsNLO.o are stored in the "include..." directories
# they are not being removed by "make clean"
# their compilation takes long so it should be done only once
#
#--------------------------------------------------------------------------------------------------



all:  $(JHUGenMELAObj) $(VegasObj) $(RockyObj) $(YetiObj) $(PDFObj) $(allObjects)
	@echo " linking"
ifeq ($(UseLHAPDF),Yes)
	echo " interfaced with LHAPDF"
else
	echo " using internal PDF sets" 
endif	
	@echo " executable file is " $(Exec)
	@echo " "
# 	$(fcomp) -o $(Exec) $(allObjects) $(RockyObj) $(YetiObj) $(IntegralObj) $(VegasObj) $(PDFObj) $(MadGraphObj)
	$(fcomp) -o $(Exec) $(allObjects) $(RockyObj) $(YetiObj) $(IntegralObj) $(VegasObj) $(PDFObj) $(JHUGenMELAObj)
# $(ObjectDir)/fastjetfortran.o $(FJLIBS) -lstdc++    add this to above line when fastjet routines are used


clean:
	rm -f ./modules/*.mod
	rm -f ./objects/*.o



./summer: summer.f90
	@echo " compiling" $<
	$(fcomp) $< -o summer


# ------------------------------------------------------------





$(ObjectDir)/mod_Misc.o: mod_Misc.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Parameters.o: mod_Parameters.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Process.o: mod_Process.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Permutations.o: mod_Permutations.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_IntegerPartition.o: mod_IntegerPartition.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_MyRecurrence.o: mod_MyRecurrence.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_MyWeylRecurrence.o: mod_MyWeylRecurrence.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Amplitudes.o: mod_Amplitudes.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/mod_SingleTopHAmps.o: mod_SingleTopHAmps.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Amplitudes_Zprime.o: mod_Amplitudes_Zprime.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/mod_Amplitudes_eeTTB.o: mod_Amplitudes_eeTTB.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

	
$(ObjectDir)/main.o: main.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTB.o: mod_CrossSection_TTB.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTBJ.o: mod_CrossSection_TTBJ.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTBP.o: mod_CrossSection_TTBP.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/mod_CrossSection_TTBP_anomcoupl.o: mod_CrossSection_TTBP_anomcoupl.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTBETmiss.o: mod_CrossSection_TTBETmiss.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_ZprimeTTB.o: mod_CrossSection_ZprimeTTB.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTBZ.o: mod_CrossSection_TTBZ.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_TTBH.o: mod_CrossSection_TTBH.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/mod_CrossSection_TH.o: mod_CrossSection_TH.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_CrossSection_eeTTB.o: mod_CrossSection_eeTTB.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_DKIntDipoles.o: mod_DKIntDipoles.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_HadrWDecay.o: mod_HadrWDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_IntDipoles.o: mod_IntDipoles.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_Integrals.o: mod_Integrals.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_JPsiFrag.o: mod_JPsiFrag.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Kinematics.o: mod_Kinematics.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Weighting.o: mod_Weighting.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@




$(ObjectDir)/mod_SixFermionProcs2.o: mod_SixFermionProcs2.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_ExoticDecay.o: mod_ExoticDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_SpinCorrel.o: mod_SpinCorrel.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_TopDecay.o: mod_TopDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(Here)/includes_DKP/DKP1L.o: $(Here)/includes_DKP/DKP1L.f
	@echo " compiling" $<
	@echo " this typically takes 10 minutes"
	$(fcomp) -80 -I$(Here)/includes_DKP -module $(ModuleDir) -c $< -o $@


$(Here)/includes_DKJ/TopCoeffsNLO.o: $(Here)/includes_DKJ/TopCoeffsNLO.f90
	@echo " compiling" $<
	@echo " this typically takes 10 minutes"
	$(fcomp) -I$(Here)/includes_DKJ -module $(ModuleDir) -c $< -o $@


$(Here)/includes_DKJ/ATopCoeffsNLO.o: $(Here)/includes_DKJ/ATopCoeffsNLO.f90
	@echo " compiling" $<
	@echo " this typically takes 10 minutes"
	$(fcomp) -I$(Here)/includes_DKJ -module $(ModuleDir) -c $< -o $@


$(ObjectDir)/mod_TTBJ_NLODK.o: mod_TTBJ_NLODK.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_TTBJ_NLODKW.o: mod_TTBJ_NLODKW.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_TTBP_NLODK.o: mod_TTBP_NLODK.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c -fixed -80 -I$(Here)/includes_DKP $< -o $@



$(ObjectDir)/mod_WDecay.o: mod_WDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(ObjectDir)/mod_ZDecay.o: mod_ZDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@




$(ObjectDir)/mod_HDecay.o: mod_HDecay.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@



$(OPPObj): $(OPPDep) $(makeDep)
	fpp @quadprec.cfg mod_NVBasis.f90  mod_NVBasis128.f90
	fpp @quadprec.cfg mod_Residues.f90 mod_Residues128.f90
	fpp @quadprec.cfg mod_UCuts.f90    mod_UCuts128.f90
	fpp @quadprec.cfg mod_Residues_new.f90  mod_Residues128_new.f90
	fpp @quadprec.cfg mod_UCuts_new.f90     mod_UCuts128_new.f90
	@echo " compiling OPP (64bit) "
	$(fcomp) -c mod_NVBasis.f90 -o $(ObjectDir)/mod_NVBasis.o
	$(fcomp) -c mod_Residues.f90 -o $(ObjectDir)/mod_Residues.o
	$(fcomp) -c mod_Residues_new.f90 -o $(ObjectDir)/mod_Residues_new.o
	$(fcomp) -c mod_UCuts.f90 -o $(ObjectDir)/mod_UCuts.o
	$(fcomp) -c mod_UCuts_new.f90 -o $(ObjectDir)/mod_UCuts_new.o
	@echo " compiling OPP (128bit) "
	$(fcomp) -r16 -c mod_NVBasis128.f90 -o $(ObjectDir)/mod_NVBasis128.o
	$(fcomp) -r16 -c mod_Residues128.f90 -o $(ObjectDir)/mod_Residues128.o
	$(fcomp) -r16 -c mod_Residues128_new.f90 -o $(ObjectDir)/mod_Residues128_new.o
	$(fcomp) -r16 -c mod_UCuts128.f90 -o $(ObjectDir)/mod_UCuts128.o
	$(fcomp) -r16 -c mod_UCuts128_new.f90 -o $(ObjectDir)/mod_UCuts128_new.o
	rm mod_NVBasis128.f90 mod_Residues128.f90 mod_UCuts128.f90 mod_UCuts128_new.f90 mod_Residues128_new.f90



$(ObjectDir)/mod_Dipoles%.o: $(DipoleDir)/mod_Dipoles%.f90 $(makeDep)
	@echo " compiling dipoles" mod_Dipoles$*".f90"
	$(fcomp) -c $(DipoleDir)/mod_Dipoles$*".f90" -o $(ObjectDir)/mod_Dipoles$*".o"


$(ObjectDir)/mod_IntDipoles%.o: $(DipoleDir)/mod_IntDipoles%.f90 $(makeDep)
	@echo " compiling dipoles" mod_IntDipoles$*".f90"
	$(fcomp) -c $(DipoleDir)/mod_IntDipoles$*".f90" -o $(ObjectDir)/mod_IntDipoles$*".o"


ifeq ($(useMPI),No)
$(VegasObj): $(VegasDir)/vegas.f $(VegasDir)/vegas_common.f $(makeDep)
	@echo " compiling" $<
	$(fcomp) -D_WriteTmpHisto=1 -c $(VegasDir)/vegas.f -o $(VegasObj)
else
$(VegasObj): $(VegasDir)/pvegas_mpi.c $(makeDep)
	@echo " compiling" $<
	$(ccomp) -D_WriteTmpHisto=1 -c $(VegasDir)/pvegas_mpi.c -o $(VegasObj)
endif
	
$(ObjectDir)/genps.o: $(PSDir)/genps.c $(makeDep)
	@echo " compiling" $<
	$(ccomp) -c $< -o $@

$(ObjectDir)/boost.o: $(PSDir)/boost.c $(makeDep)
	@echo " compiling" $<
	$(ccomp) -c $< -o $@

$(YetiObj): $(PSDir)/yeti.f $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $(YetiObj)
	
$(ObjectDir)/mstwpdf.o: $(PDFDir)/mstwpdf.f $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@
	
$(ObjectDir)/mrst2001lo.o: $(PDFDir)/mrst2001lo.f $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/cteq2mrst.o: $(PDFDir)/cteq2mrst.f $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/Cteq66Pdf.o: $(PDFDir)/Cteq66Pdf.f $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/CT10Pdf.o: $(PDFDir)/CT10Pdf.f $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(ObjectDir)/NNPDFDriver.o: $(PDFDir)/NNPDFDriver.f $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(JHUGenMELAObj): $(JHUGenMELADep) $(JHUGenMELADir)/variables.F90 $(JHUGenMELADir)/includeVars.F90
	@echo " compiling JHUGenMELA"
	$(fcomp) -c $(JHUGenMELADep) -o $@




# supresses command calls
.SILENT:
