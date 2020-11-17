# ubsmear

A dependency-free header-only C++ library to forward-fold cross-section predictions so they can be compared to MicroBooNE data.

## Installation

```
# Clone this repo
git clone https://github.com/a-d-smith/ubsmear.git

# Simply add ubsmear/inc as an include path, and #include "ubsmear.h"
# For a simple example, see ubsmear/example/Makefile.
```

## How to compile the example

```
# Clone this repo
git clone https://github.com/a-d-smith/ubsmear.git

# Go to the ubsmear directory
cd ubsmear

# Set the UBSMEAR_DIR environment variable to point to your clone of the repo
export UBSMEAR_DIR=`pwd`

# Go to the example directory, and compile the test code
cd example
make

# Run the example
./test
```
