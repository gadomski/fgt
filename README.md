# fgt

Fast Gauss transforms.

The Gauss transform is a common operation that computes the per-point similarity between two data sets:

![The Gauss transform](img/gauss-transform.png)

This a C++ library for computing the Gauss transform using the direct method as well as a few shortcuts.
This code lives on [Github](https://github.com/gadomski/fgt), has [some Doxygen documentation](http://gadomski.github.io/fgt), and is tested with [Travis](https://travis-ci.org/gadomski/fgt) and [AppVeyor](https://ci.appveyor.com/project/gadomski/fgt/branch/master).

![Travis Build Status](https://travis-ci.org/gadomski/fgt.svg?branch=master)
![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/7t6ojbj2lj052wx2/branch/master?svg=true)


## Usage

There is one C++ header file, `fgt.hpp`, which has everything you need.
Include that file and you're off to the races:

```cpp
#include <fgt.hpp>

void my_great_function(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y) {
    double bandwidth = 0.3;
    Eigen::VectorXd gauss_transform = fgt::direct(x, y, bandwidth);
}
```

The library provides a few different ways to calculate the Gauss transform:

- `fgt::direct` calculates the exact Gauss transform, and is the most accurate and the slowest option.
- `fgt::direct_tree` tries to do less work by only considering "close" points, where "close" is defined by the bandwidth.
  The direct tree method works best for very small bandwidths.
- `fgt::ifgt` uses the [Improved Fast Gauss Transform (pdf)](http://www.umiacs.umd.edu/~yangcj/papers/siam_fgt_v11.pdf) to speed up the calculation.
  IFGT is fast for large bandwidths but can break down for smaller bandwidths.

There's also a class-based interface:

```cpp
#include <fgt.hpp>

void my_great_function(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y) {
    double bandwidth = 0.3;
    fgt::Direct direct(x, bandwidth);
    Eigen::VectorXd result = direct.compute(y);
}
```

This lets you break up your transform into a pre-compute and a compute step, which can save you some cycles if you're re-using the same source dataset in a more advanced transform (e.g. direct_tree or ifgt).

There is some benchmarking code available in the [bench](bench/) directory, which you can use to try to get a sense of the performance of the various modes.
We found a crossover point at bandwidths of a bit more than 0.2 during local testing on a Mac laptop; YMMV.

![Benchmarks conducted on a random dataset on my personal Mac laptop](img/2016-03-01-clang-700.1.81.x86_64-apple-darwin15.3.0.png)

## Installation

**fgt** has no runtime dependencies, and only depends on [CMake](https://cmake.org/) and [Eigen](http://eigen.tuxfamily.org/) for building.

### Homebrew

If you're on a Mac and you use [Homebrew](http://brew.sh), use [my tap](https://github.com/gadomski/homebrew-gadomski) to install:

```sh
brew tap gadomski/gadomski
brew install fgt
```

### From source

To build **fgt** from source, clone the repository and execute the traditional CMake build incantation:

```sh
git clone https://github.com/gadomski/fgt
mkdir -p fgt/build
cd fgt/build
cmake ..
make
```

**fgt** doesn't make any assumptions about whether you do or do not want shared libraries, so if you have a preference be sure to set `BUILD_SHARED_LIBS`.

### Eigen ordering

Eigen, by default, stores matrices in column-major order, but **fgt** works with [row-major](https://en.wikipedia.org/wiki/Row-major_order) matrices.
If you want to avoid extra copies, pass in row-major matrices to **fgt** functions.
You can use the `fgt::Matrix` typedef to help:

```cpp
fgt::Matrix my_matrix(1000, 3); // creates an uninitialized 1000x3 row-major matrix of doubles
```


### OpenMP

**fgt** comes with built-in OpenMP parallelization, which can lead to some significant speedups for large data sets.
To enable OpenMP, make sure you're using an OpenMP-aware compiler (on OSX, you can get OpenMP clang via Homebrew: `brew install clang-omp`) and set the CMake variable `WITH_OPENMP` to ON, e.g.:

```sh
CC=clang-omp CXX=clang-omp++ cmake .. -DWITH_OPENMP=ON
make
```

This will build an OpenMP-enabled **fgt** library.

### Tests

**fgt** comes with a unit-test suite.
To run, simply execute `make test`.

### Using in a downstream project

When you install **fgt** on your system, you will also install a few CMake configuration files that make it easy to integrate this project into your downstream work.
If **fgt** is installed to a traditional location (e.g. `/usr/local`), finding **fgt** might be as simple as including the following in your `CMakeLists.txt`:

```cmake
find_package(Fgt REQUIRED)
target_link_libraries(my-sweet-target
    PUBLIC
    Fgt::Library-C++
    )
```

The provided target, `Fgt::Library-C++`, should have its `INTERFACE_INCLUDE_DIRECTORIES` and other useful properties configured, so you shouldn't have to muck with anything else.

If you've installed **fgt** to a non-standard location, you may need to use `Fgt_DIR` to find its CMake configuration files, or use `CMAKE_PREFIX_PATH` to point CMake at your non-standard install tree.

## Versioning

We follow [semantic versioning](http://semver.org/).
Versions have annotated tags following a `vMAJOR.MINOR.PATCH` naming convention.
While we'll do our best to increment the `MINOR` version with all breaking changes, we can't guarantee anything until `MAJOR` hits 1 (per usual).

## Contributing

As always, we welcome bug reports, features requests, and particularly pull requests via [Github](https://github.com/gadomski/fgt).

This library was developed by [Pete Gadomski](https://github.com/gadomski), and it was inspired by [the IFGT source code](http://www.umiacs.umd.edu/labs/cvl/pirl/vikas/Software/IFGT/IFGT_code.htm) and [figtree](https://github.com/vmorariu/figtree).

## License

LGPL v2.1.
See LICENSE.txt for the complete text.
