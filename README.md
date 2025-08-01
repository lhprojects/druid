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

If LCIO can not be found in the env, we will download and build it automatically.

**Note:**
If you met issue about `fdopen`, try commenting out line `define fdopen(fd,mode) NULL` in 
`build/_deps/lcio/src/LCIO_ext-build/_deps/sio_extern-src/zlib/zutil.h`.
This solves problem, when your system has `fdopen`, but it thinks your system has no one.

### cmake
`cmake` version 3.30 works well, but version 4 is too high.

## Build and Install Druid

### Install Druid with CMake
#### Prepare
setup ROOT, so that the build system can see it.
Normally, you don't have to do anything, if `root` command works.
#### Builds

   ```
   mkdir build && cd build
   cmake ..
   make && make install
   ```

   CMAKE will detect the platform, search the path of dependencies, and generate setup.sh file automatically.
   After compiling, the setup.sh of Druid and executable file will be generated in <path-druid>/bin.
   
   
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
