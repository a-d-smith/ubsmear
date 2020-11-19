# Unit tests

This directory contains a suite of unit tests for the main code base.

## How to run the tests

```
# Make sure the UBSMEAR_DIR environment variable to point to your clone of the ubsmear repo
export UBSMEAR_DIR=/path/to/ubsmear/

# Compile and run all of the tests
cd $UBSMEAR_DIR/test/
./runTests
```

## How to write a new test

```
# Think of a name for your test, here we will use "example"
# Make a directory for the test and prefix it with "test_"
cd $UBSMEAR_DIR/test/
mkdir test_example

# Make a new file called test.cc, which has an `int main()` method (or copy the basic test and edit that)
cp test_basic/test.cc test_example/.

# Write (or copy) a Makefile
cp test_basic/Makefile test_example/.

# If you used the above naming structure, it should be picked up when you execute ./runTests
```

## How to run an individual test
```
# Go to the directory of the test you want to run
cd $UBSMEAR_DIR/test/test_basic

# Compile the test (you might also want to run `make clean` if you have changed something in ubsmear but not the test itself)
make

# Run the test and print the exit code
./test; echo "Exit code: $?"

# If the exit code is zero, then the test passed!
```
