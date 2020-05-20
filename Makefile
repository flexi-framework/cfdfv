### CFDFV Makefile

include config.mk

.PHONY: all cfdfv shared clean veryclean cleanshare release

all: shared cfdfv
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: ALL EXECUTABLES GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

cfdfv:
	cd src/ && touch deplist.mk && $(MAKE) all

shared:
	cd share/ && $(MAKE) all

clean:
	cd src/ && $(MAKE) clean

veryclean:
	cd src/ && $(MAKE) veryclean

cleanshare:
	cd share/ && $(MAKE) clean

release:
	mkdir CFDFV
	cp -r bin CFDFV
	cp -r Calc CFDFV
	cp -r lib CFDFV
	cp -r share CFDFV
	cp -r src CFDFV
	cp -r utils CFDFV
	cp Makefile config.mk CFDFV
	/bin/sh utils/relcleanup.sh CFDFV
	tar cjf CFDFV.tar.bz2 CFDFV
