all:
	$(MAKE) -C build all

debug:
	$(MAKE) -C build debug   

classic:
	$(MAKE) -C build classic

debug_classic:
	$(MAKE) -C build debug_classic

gfortran:
	$(MAKE) -C build gfort

debug_gfortran:
	$(MAKE) -C build debug_gfort

clean:
	$(MAKE) -C build clean

