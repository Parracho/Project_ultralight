see which version of gsl: $gsl-config --version

to create the "build" folder: Cuba-4.2/codefolder/$ mkdir build

in the build folder:	$ cmake ..

	files created in build: CMakeCache.txt  CMakeFiles  cmake_install.cmake  DM_freeze_in  Makefile

			$ make

			$ cmake -DCMAKE_INSTALL_PREFIX=.. ..

			$ make

			$ make install

in the bin folder:

			$./DMFI	