
default: all

all: qe2pert perturbo

qe2pert: pw ph
	cd qe2pert-src ; $(MAKE)
	(mkdir -p bin)
	(cd bin; ln -sf ../qe2pert-src/qe2pert.x .)

perturbo: pw ph
	cd pert-src ; $(MAKE)
	(mkdir -p bin)
	(cd bin; ln -sf ../pert-src/perturbo.x .)

tdep:
	cd external/tdep && make

pw: 
	cd .. ; $(MAKE) pw

ph: 
	cd .. ; $(MAKE) ph

clean:
	cd pert-src ; $(MAKE) clean
	cd qe2pert-src ; $(MAKE) clean
	cd external/tdep && make clean

.PHONY: default all clean pw ph qe2pert perturbo
