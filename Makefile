all:
	$(MAKE) -C build all

debug:
	$(MAKE) -C build debug   

gfortran:
	$(MAKE) -C build gfort

debug_gfortran:
	$(MAKE) -C build debug_gfort

clean:
	$(MAKE) -C build clean

