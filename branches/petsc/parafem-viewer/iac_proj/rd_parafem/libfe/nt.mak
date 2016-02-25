# Build the felib.dll library

all:
	cd src && $(MAKE) -f nt.mak

clean:
	cd src && $(MAKE) -f nt.mak clean

