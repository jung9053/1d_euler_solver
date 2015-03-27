include ./Makefile.in

MYNAME = eu1
MYNAME1 = EU1
VER = 10
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
DATE = $(shell date +%Y%m%d.%HH)
BACKUPNAME = $(MYNAME1)_$(VER)_$(DATE).tar

#LFLAGS = ./tecio.a -lstdc++ -L. -linvmat
LFLAGS =

FFLAGS = $(FOPTFLAGS)
CFLAGS = $(COPTFLAGS) -I. $(INCDIR)

#----------------------------------------------------------------------
#main routine
FSRC0 = drive_euler_1d.f
#----------------------------------------------------------------------
#when the source should be passwd through C preprocessor
FSRC1 =
#----------------------------------------------------------------------
#when the source don't need to be passwd through C preprocessor
FSRC2 = module_1d.f\
	\
	euler_1d.f	rinput_1d.f	set_dimen_1d.f\
	allocate_1d.f	make_grid_1d.f	grid_metric_1d.f\
	initialize_1d.f	\
	qface_1d.f	calc_roe_flux_1d.f	conpq_1d.f\
	calc_inv_flux_1d.f	uctime_1d.f\
	main_iter_1d.f	calc_res_1d.f\
	euler_fwd_1d.f\
	wrt_out_1d.f  muscl_recon.f  weno5_recon.f crweno5_recon.f tri.f\
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#when the source should be passwd through C preprocessor
#CSRC1=
#----------------------------------------------------------------------
#when the source don't need to be passwd through C preprocessor
CSRC2=
#----------------------------------------------------------------------

FSRC = $(FSRC0) $(FSRC1) $(FSRC2)

CSRC = $(CSRC1) $(CSRC2)

SRC  = $(FSRC) $(CSRC)

ifneq ($(all),1)
OBJS = $(CSRC:.c=.o) $(FSRC1:.F=.o) $(FSRC2:.f=.o) $(FSRC0:.f=.o)
else
OBJS = $(CSRC:.c=.o) $(FSRC1:.F=.o) $(FSRC2:.f=.o)
endif

BACKUPFILE = $(FSRC) euler.inp Makefile Makefile.in ffscan

.F.o :
	$(FC) $(FFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<
.c.o :
	$(CC) $(CFLAGS) -c $<

default : $(where)

mk_exe: $(OBJS)
	$(FC) -o $(EXCNAME) $(OBJS) $(LFLAGS)

mk_lib: $(OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBNAME)

clean:
	rm -f *~
	rm -f *.o 
	rm -f *.mod
	rm -f *__genmod.f90
	rm -f $(EXCNAME)
	rm -f $(LIBNAME)
	rm -f $(SRCNAME)
	
new : clean default

src :
	cat $(SRC) > $(SRCNAME)
tag :
	ctags $(SRC)
backup :
	tar cvf $(BACKUPNAME) $(BACKUPFILE)

scan :
	./ffscan $(FSRC)

checkin:
	@for file in *.[c,h]; \
	do \
	ci -u -m'Maintance' $$file;\
	done 

checkin2:
	@for file in *.[c,h]; \
	do \
	ci $$file;\
	rcs -U $$file;\
	co $$file;\
	done
