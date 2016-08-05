CHARMC= /projects/ATPESC2016/anl_charm_tutorial/charm/bin/charmc $(OPTS)

OBJS = keyvalue.o

all: keyvalue

keyvalue: $(OBJS)
	$(CHARMC) -language charm++ -o keyvalue $(OBJS)

keyvalue.decl.h: keyvalue.ci
	$(CHARMC)  keyvalue.ci

clean:
	rm -f *.decl.h *.def.h conv-host *.o keyvalue charmrun

keyvalue.o: keyvalue.C keyvalue.decl.h
	$(CHARMC) -c keyvalue.C

test: all
	./charmrun +p8 ++local keyvalue 
