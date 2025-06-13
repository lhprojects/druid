# Install LCIO and Druid on MacOS

## Software version

- MacOS: 

  系统版本：	macOS 12.3.1 (21E258)
  
  内核版本：	Darwin 21.4.0

- ROOT:

  ```bash
  ➜  README git:(dev) ✗ root
     ------------------------------------------------------------------
    | Welcome to ROOT 6.24/06                        https://root.cern |
    | (c) 1995-2021, The ROOT Team; conception: R. Brun, F. Rademakers |
    | Built for macosx64 on Sep 02 2021, 14:20:23                      |
    | From tags/v6-24-06@v6-24-06                                      |
    | With Apple clang version 12.0.5 (clang-1205.0.22.9)              |
    | Try '.help', '.demo', '.license', '.credits', '.quit'/'.q'       |
     ------------------------------------------------------------------
  ```

  v6.26.02 does not work on macOs Monterey.

  v6.24.06 built with conda, cvmfs and pre-compiled binary distribution are all tested.

- LCIO: v02-17

- Druid: Druid_2.4

LCIO and Druid are included in this repository.

## Dependency

### ROOT

See [ROOT](https://root.cern/install/)

> Some versions of pre-compiled binary do not include all components needed,
so it may be the easiest way to [install root from conda](https://root.cern/install/#conda).

### Install LCIO

In IHEP cluster, the LCIO and ROOT are built in [LCG release](https://lcgdocs.web.cern.ch/lcgdocs/lcgreleases/introduction/),
e.g. users can go into LCG100 by `source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc11-opt/setup.sh`.

If users do not have LCIO in their device, they can install it from the official webpage, 
[LCIO (desy.de)](https://ilcsoft.desy.de/portal/software_packages/lcio/), or use the attached version.

1. Install java

   ```bash
   brew install java
   ```

2. Install LCIO

   ```bash
   # For instance, assuming this repository is installed as $HOME/opt/druid.
   cd ~/opt
   git clone https://code.ihep.ac.cn/cheyuzhi/druid.git
   cd druid/LCIO
   
   mkdir build
   cd build
   cmake ..
   make -j 8 install
   ```

> Note: Under Fedora, zlib should be installed by yum as a requirement:
> ```bash
> sudo yum install zlib zlib-devel
> ```

## Install Druid
### Install Druid with GNUMakefile

   before that, activate your conda environment of root, or source `thisroot.sh` in `<root-location>/bin`.
   
   And edit variable `${LCIO}` and `${DRUIDDIR}` in druid/Druid_2.4/setup.sh:
   
   ```bash
   # Except following variables, one should source thisroot.sh
   unset LCIO
   export LCIO=$HOME/opt/druid/LCIO
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LCIO/lib
   
   export DRUIDDIR=$HOME/opt/druid/Druid-2.4
   
   export PATH=$PATH:$DRUIDDIR:$LCIO/bin:$DRUIDDIR/bin 
   
   ```

   > Under Fedora, or maybe other Linux distributions, the LCIO library is installed in ${LCIO}/lib64. Therefore
   the LCIO_LIBRARY_PATH should be changed to
   > ```bash
   >    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LCIO/lib64
   > ```
   
   Then compile Druid:
   
   ```bash
   cd druid/Druid_2.4
   . ./setup.sh
   
   cd src
   make -j 8
   ```

### Install Druid with CMake
   
    We provide a CMakeLists.txt file, which showed better stability under newer ROOT & LCIO softwares and multiple platforms.
    Now this approach has been tested under IHEP cluster (with CentOS7 + LCG100) and MacOS 12.6.6.

   To compile Druid with CMake, one only needs:

   1. ensure that cmake works,
   2. get into ROOT and LCIO surrounding
   
   and then,

   ```
   mkdir build && cd build
   cmake ..
   make && make install
   ```

   CMAKE will detect the platform, search the path of dependencies, and generate setup.sh file automatically.
   After compiling, the setup.sh of Druid and executable file will be generated in <path-druid>/bin.

   NOTE: There are two `setup.sh` files, the one mentioned in GNUmakefile approach should be modified properly by users, 
   while the another one is generated from <path-druid>/setup.sh.in by CMAKE and located in <path-druid>/bin. We keep them in 
   the same while, so that the two approaches work in paralell.
   
4. Enjoy :)

	```bash
	Druid <an_slcio_sample>
	```