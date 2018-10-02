#
# object files
# ------------
#
OBJ1 = drop_3dw.o trgl6_hsph_octa.o crvm_3d.o crvm_3d_interp.o
OBJ2 = sgf_3d_w.o sgf_3d_fs.o
OBJ3 = dfel_3d.o dfel_3d_interp.o 
OBJ39= abc.o elm_geom.o printel.o inclination.o taylor.o xy_slice.o xz_slice.o interp_p.o
OBJ4 = sdlp_3d.o sdlp_3d_interp.o sdlp_3d_integral.o
OBJ41= sslp_3d.o sslp_3d_integral.o sslp_3d_integral_sing.o sslp_3d_interp.o
OBJ5 = vel.o deflation.o
OBJ6 = cramer_33.o gauss_trgl.o gauss_disk.o gauss_leg.o
OBJ7 = sgrad_3d.o sgrad_3d_interp.o
OBJA =  $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ41)
OBJB =  $(OBJ39)
OBJ  =  $(OBJA) $(OBJB)
# 
# link
# ----
#
drop_3dw: $(OBJ)
	gfortran -o drop_3dw $(OBJ) -llapack -lblas
#
# compile
# -------
#
drop_3dw.o: drop_3dw.f
	gfortran -c drop_3dw.f
trgl6_hsph_octa.o: trgl6_hsph_octa.f 
	gfortran -c trgl6_hsph_octa.f 

sdlp_3d.o: sdlp_3d.f 
	gfortran -c sdlp_3d.f   
sdlp_3d_interp.o: sdlp_3d_interp.f 
	gfortran -c sdlp_3d_interp.f   
sdlp_3d_integral.o: sdlp_3d_integral.f 
	gfortran -c sdlp_3d_integral.f  
gauss_disk.o: gauss_disk.f 
	gfortran -c gauss_disk.f 
elm_geom.o: elm_geom.f
	gfortran -c elm_geom.f   
abc.o: abc.f
	gfortran -c abc.f  
printel.o: printel.f
	gfortran -c printel.f  
interp_p.o: interp_p.f
	gfortran -c interp_p.f 
taylor.o: taylor.f
	gfortran -c taylor.f 
inclination.o: inclination.f
	gfortran -c inclination.f 
#
#
# all below borrowed from drop_3d
# -------------------------------
#
#
 vel.o: vel.f 
	gfortran -c  vel.f 
deflation.o: deflation.f
	gfortran -c deflation.f 
dfel_3d.o: dfel_3d.f 
	gfortran -c dfel_3d.f   
dfel_3d_interp.o: dfel_3d_interp.f 
	gfortran -c dfel_3d_interp.f  
sslp_3d.o: sslp_3d.f
	gfortran -c sslp_3d.f  
sslp_3d_interp.o: sslp_3d_interp.f
	gfortran -c sslp_3d_interp.f  
sslp_3d_integral.o: sslp_3d_integral.f
	gfortran -c sslp_3d_integral.f 
sslp_3d_integral_sing.o: sslp_3d_integral_sing.f
	gfortran -c sslp_3d_integral_sing.f  
crvm_3d.o: crvm_3d.f
	gfortran -c crvm_3d.f
crvm_3d_interp.o: crvm_3d_interp.f
	gfortran -c crvm_3d_interp.f 
xy_slice.o: xy_slice.f
	gfortran -c xy_slice.f 
xz_slice.o: xz_slice.f
	gfortran -c xz_slice.f
sgrad_3d.o: sgrad_3d.f
	gfortran -c sgrad_3d.f
sgrad_3d_interp.o: sgrad_3d_interp.f
	gfortran -c sgrad_3d_interp.f
cramer_33.o: cramer_33.f
	gfortran -c cramer_33.f 
sgf_3d_w.o: sgf_3d_w.f 
	gfortran -c sgf_3d_w.f 
sgf_3d_fs.o: sgf_3d_fs.f
	gfortran -c sgf_3d_fs.f 
gauss_leg.o: gauss_leg.f 
	gfortran -c gauss_leg.f 
gauss_trgl.o: gauss_trgl.f 
	gfortran -c gauss_trgl.f
#
# all
# ---
#
all:
	make drop_3dw
#
# clean
# -----
#
clean:
	rm -f core $(OBJ) drop_3dw
	rm -f drop_3dw.xyz drop_3dw.xy drop_3dw.diag drop_3dw.net
#
# tar
# ---
#
tarit:
	tar -cvf drop_3dw.tar makefile *.f *.dat
