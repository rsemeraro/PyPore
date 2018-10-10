[:arrow_left:](https://github.com/rsemeraro/PyPore/blob/master/README.md)

## Installation of OpenMPI on MAC OS X

Message Passing Interface (MPI) is an application interface that allows processes to communicate by sending and receiving messages. In essence, it facilitates communication when multiple computers/processors/cores are performing distributed calculations. For people familiar with computer clusters, MPI is the interface that helps the compute nodes talk to each other in parallel computations. I won't provide detailed information about MPI and programming practices for using MPI here, but rather the purpose of this page is to help you install one of the freely available MPI implementations ( **Open MPI** ) on a computer running Mac OS X.

You can download the last version from it's official page:

```
http://www.open-mpi.org/software/
```
**Prerequisites:**

To install **Open MPI** as described in this HOWTO, you'll need the following:
- A Mac running Mac OS 10.6 (Snow Leopard), 10.7 (Lion), 10.8 (Mountain Lion) or above
- Internet access
- Fortran and or C/C++ compilers

### Manual install

The installation of **Open MPI** on Mac is surprisingly easy. If have not already installed the GNU compilers, check out how to install the GNU compilers on Mac OS X.

1. Copy the archive file to a temporary location to compile it. Open **Terminal.app** and change directories to the new location of the **Open MPI** archive. If you downloaded the source code as a .tar.gz file, you can untar/zip it by typing:
    ```
    tar zxvf openmpi-1.6.5.tar.gz
    ```
    which will create a directory called `openmpi-1.6.5/`. Change directories into the new directory.
2. You can now run the configuration script. If you only have one set of compilers installed, you can run the configuration script by typing:
    ```
    ./configure --prefix=/usr/local
    ```
    If you would like to install Open MPI someplace other than `/usr/local`, you can change the installation directory listed after the prefix flag.
    If you have multiple compilers installed, you can specify which compilers you would like to use as follows:
    ```
    ./configure CC=icc CXX=icpc F77=ifort FC=ifort --prefix=/usr/local
    ```
    where you specify the C (CC), C++ (CXX), Fortran 77 (F77) and Fortran 90 and later (FC) compilers with the listed variables.
3. Assuming the configuration script runs without issue, you can compile the code by typing:
    ```
    make all
    ```
    which will compile the Open MPI libraries/binaries and configure the wrappers for using the specified compilers. This should take a bit...

4. Again, assuming the code compiles without issue, you can install **Open MPI** by typing:
    ```
    sudo make install
    ```
    Beware that using sudo can do major damage to your computer if you aren't careful. You can now feel free to delete the temporary Open MPI build directory (e.g., openmpi-1.6.5/) since the code has been installed. If you think you might want to rebuild with different compilers or change the configuration, you may want to keep this directory around.

So now you should have a functional installation of **Open MPI** on your Mac. You should be able to compile code that uses MPI by using the **Open MPI** compiler wrappers (mpicc, mpicxx, mpif77, mpif90) and run MPI-enabled programs with mpiexec. If you try to use the new **Open MPI** executables and they are not found, it may be that /usr/local/bin (or wherever you specified with the --prefix flag in the configure stage) is not in your $PATH environment variable. You can confirm this by typing:
```
echo $PATH
```

### Some common errors

If you get an error like the next, during PyPore installation:
```
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
warning: build_clib: command '/usr/local/bin/mpicc' failed with exit status 1
warning: build_clib: building optional library "mpe" failed
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
warning: build_clib: command '/usr/local/bin/mpicc' failed with exit status 1
warning: build_clib: building optional library "vt" failed
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
warning: build_clib: command '/usr/local/bin/mpicc' failed with exit status 1
warning: build_clib: building optional library "vt-mpi" failed
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
warning: build_clib: command '/usr/local/bin/mpicc' failed with exit status 1
warning: build_clib: building optional library "vt-hyb" failed
ld: warning: The i386 architecture is deprecated for macOS (remove from the Xcode build setting: ARCHS)
gcc: error: unrecognized command line option '-Wshorten-64-to-32'
error: Setup script exited with error: Cannot compile MPI programs. Check your configuration!!!
```
It's because gcc is not installed in /usr/local/bin by xcode. So CC=clang is the simplest way to answer it. Type:
```
CC=/usr/bin/clang CFLAGS="-O" python setup.py install
```
If you get some error related to hdf5 library missing, likew following:
```
./h5py/api_compat.h:27:10: fatal error: 'hdf5.h' file not found
```
Probably you have to install the hdf5 library o your mac. You can obtain it by using mac brew:
```
brew install homebrew/core/hdf5
```
