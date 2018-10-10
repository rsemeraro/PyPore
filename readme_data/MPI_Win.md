[:arrow_left:](https://github.com/rsemeraro/PyPore/blob/master/README.md)
---
## Installation of OpenMPI on Windows 7/8/10
Message Passing Interface (MPI) is an application interface that allows processes to communicate by sending and receiving messages. In essence, it facilitates communication when multiple computers/processors/cores are performing distributed calculations. For people familiar with computer clusters, MPI is the interface that helps the compute nodes talk to each other in parallel computations. I won't provide detailed information about MPI and programming practices for using MPI here, but rather the purpose of this page is to help you install one of the freely available MPI implementations ( **Open MPI** ) on a computer running Windows 7/8/10.

### Manual install
Microsoft has an implementation of MPI, called as MS MPI. Download MS MPI packages from the official Microsoft website:
```
https://www.microsoft.com/en-us/download/details.aspx?id=56511
```
1. Typically, you have to download MSMPI SDK which has set of libraries to compile and linkup MPI applications on Windows and MSMPI executables to launch and run MPI applications:
    - msmpisdk.msi
    - msmpisetup.exe


2. Install both of them. Just double click and accept license agreement. Installation is pretty simple.

3. Type `mpiexec` in a terminal shell to check MPI installation 
