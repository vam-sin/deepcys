#
#	makefile for espot
#
#	april 15, 1991, e.l. mehler
#       revised May 5, 1997
#	env.f added Oct 13, 1997
#	dept. of physiology/biophysics, Mt. Sinai School of Med.
#	New York
#


FFO = espmain.o filopn.o coorin.o  cntrl.o  resfsa.o env.o rdchg.o rdsenv.o rlen.o gepol.o 
COORO = micro_debug_$(HOSTTYPE)
FC = fort77
# FC = f77
FLAGS =	-O2 -static #-fsyntax-only
.SUFFIXES: .inc

micro_$(HOSTTYPE) : $(FFO) 
	$(FC)  $(FLAGS) -o $(COORO) $(FFO) 

espmain.o : espmain.f
	$(FC) -c $(FLAGS) espmain.f -o espmain.o

cntrl.o     : cntrl.f 
	$(FC) -c $(FLAGS) cntrl.f -o cntrl.o

filopn.o : filopn.f
	$(FC) -c $(FLAGS) filopn.f -o filopn.o

rlen.o : rlen.f
	$(FC) -c $(FLAGS) rlen.f -o rlen.o

coorin.o : coorin.f
	$(FC) -c $(FLAGS) coorin.f -o coorin.o
gepol.o : gepol.f
	$(FC) -c $(FLAGS) gepol.f -o gepol.o



resfsa.o : resfsa.f
	$(FC) -c $(FLAGS) resfsa.f -o resfsa.o

rdchg.o : rdchg.f
	$(FC) -c $(FLAGS) rdchg.f -o rdchg.o

rdsenv.o : rdsenv.f
	$(FC) -c $(FLAGS) rdsenv.f -o rdsenv.o

env.o :    env.f
	$(FC) -c $(FLAGS)    env.f -o    env.o

clean :
	rm $(COORO) $(FFO)

print :
	lp $(FF)

print :
	lp $(FCOOR)
