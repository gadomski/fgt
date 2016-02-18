build:
	mkdir build; cd build; \
		cmake .. -DWITH_TESTS=ON -DWITH_OPENMP=$(FGT_WITH_OPENMP) -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON; \
		make

test:
	build/bin/fgt-test
	build/bin/fgt-trial

install:
	cd build; sudo make install


.PHONY: build test install
