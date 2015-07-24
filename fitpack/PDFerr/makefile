FC = gfortran

part =  ../Apps/pdf10.o \
	../theory/altpar10.o ../theory/xpij2p.o \
	../theory/xpns2p.o \
	../util/polint.o ../util/alpha_cwz.o ../util/libmathlib.a \
	../util/io.o ../util/trmstr.o ../util/linalg.o ../util/polint4.o \

PDFerr : PDFerr.o $(part)
	$(FC) -g -o ./$@ PDFerr.o $(part)
