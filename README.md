# Build and Install Druid On Unix-Like System

## Dependency

### ROOT

See [ROOT](https://root.cern/install/)

> Some versions of pre-compiled binary do not include all components needed,
so it may be the easiest way to [install root from conda](https://root.cern/install/#conda).

ROOT from `source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-01-28` is tested.
```bash
$ root --version
ROOT Version: 6.32.04
Built for linuxx8664gcc on Oct 28 2024, 12:41:01
From heads/master@tags/v6-32-04
```

### LCIO

LCIO from `source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-01-28` is tested.


## Build and Install Druid

### Install Druid with CMake (Suggested)

   ```
   mkdir build && cd build
   cmake ..
   make && make install
   ```

   CMAKE will detect the platform, search the path of dependencies, and generate setup.sh file automatically.
   After compiling, the setup.sh of Druid and executable file will be generated in <path-druid>/bin.
   
### Install Druid with GNUMakefile

⚠️ **Note:** This is outdated.

   before that, activate your conda environment of root, or source `thisroot.sh` in `<root-location>/bin`.
   
   And edit variable `${LCIO}` and `${DRUIDDIR}` in druid/Druid/setup.sh:
   
   ```bash
   # Except following variables, one should source thisroot.sh
   unset LCIO
   export LCIO=$HOME/opt/druid/LCIO
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LCIO/lib
   
   export DRUIDDIR=$HOME/opt/druid/Druid
   
   export PATH=$PATH:$DRUIDDIR:$LCIO/bin:$DRUIDDIR/bin 
   
   ```

   > Under Fedora, or maybe other Linux distributions, the LCIO library is installed in ${LCIO}/lib64. Therefore
   the LCIO_LIBRARY_PATH should be changed to
   > ```bash
   >    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LCIO/lib64
   > ```
   
   Then compile Druid:
   
   ```bash
   cd druid/Druid
   . ./setup.sh
   
   cd src
   make -j 8
   ```


   
# Usage
A single `Druid` without arguments to see usage.
```bash
Druid
```
An example to use
```bash
Druid <an_slcio_sample> <evtNum>
```

# Troubleshotting
If your opengl crashed, try using software rendering
```bash
export LIBGL_ALWAYS_SOFTWARE=1
Druid ...
```
