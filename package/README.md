# MATLAB Interface Packaging Functions

Make sure you have `cmake` (3.15+) on your path. No submodules are needed — OSQP v1.0.0
is fetched automatically via CMake FetchContent during the build.

Run `package_osqp.m` from within MATLAB:

```matlab
cd package
package_osqp
```

This will compile the interface and package it as a `osqp-matlab-<platform>64.tar.gz` file.

From the command line:

```
/path/to/matlab -nodisplay -nosplash -nodesktop -r "cd package; package_osqp(); exit;"
```

You can pass a version string to override the auto-detected version:

```
/path/to/matlab -nodisplay -nosplash -nodesktop -r "cd package; package_osqp('1.0.0'); exit;"
```

Once the `.tar.gz` files for each platform have been generated, upload them to the
appropriate GitHub release as assets.
