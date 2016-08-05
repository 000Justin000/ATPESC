OPTS=-O3
CHARMC= /projects/ATPESC2016/anl_charm_tutorial/charm/bin/charmc $(OPTS)

all: leanmd

leanmd: Main.o Cell.o Compute.o leanmd.decl.h
	$(CHARMC) $(OPTS) -module DummyLB \
	-language charm++ -tracemode projections -o leanmd Main.o Cell.o Compute.o

Main.o: Main.C Main.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o $@ $<

Cell.o: Cell.C Cell.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o $@ $<

leanmd.decl.h:	leanmd.ci
	$(CHARMC) leanmd.ci

Compute.o: Compute.C Compute.h leanmd.decl.h defs.h physics.h
	$(CHARMC) $(OPTS) -o $@ $<

test: leanmd
	$(call run, +p4 ./leanmd 3 3 3 20 )

clean:
	rm -f *.decl.h *.def.h *.o leanmd leanmd.prj charmrun
