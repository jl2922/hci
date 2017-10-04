#!/bin/bash
# Install dependencies and run tests.

test -n $CC && unset CC
test -n $CXX && unset CXX

# Install or Load OpenMPI.
if [ -f "$TOOLS_DIR/openmpi/bin/mpic++" ] && [ -f "$TOOLS_DIR/openmpi/bin/mpic++" ]; then
	echo "Found cached OpenMPI"
else
	echo "Downloading OpenMPI Source"
  mkdir -p downloads
  cd downloads
	wget -O openmpi-3.0.0.tar.bz2 https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.bz2
	tar xjf openmpi-3.0.0.tar.bz2
	echo "Configuring and building OpenMPI"
	cd openmpi-3.0.0
  mkdir -p $TOOLS_DIR/openmpi
	./configure --prefix=$TOOLS_DIR/openmpi CC=$C_COMPILER CXX=$CXX_COMPILER
	make -j 8
	make install
	echo "Completed"
	echo
	cd ../../
fi
export PATH=$TOOLS_DIR/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$TOOLS_DIR/openmpi/lib:$LD_LIBRARY_PATH

# Install or Load Protocol Buffers.
if [ -f "$TOOLS_DIR/protobuf/bin/protoc" ]; then
	echo "Found cached Protocol Buffers"
else
	echo "Downloading Protocol Buffers"
  mkdir -p downloads
  cd downloads
	wget -O protobuf-cpp-3.4.1.tar.gz https://github.com/google/protobuf/releases/download/v3.4.1/protobuf-cpp-3.4.1.tar.gz
	tar xzf protobuf-cpp-3.4.1.tar.gz
	echo "Configuring and building Protocol Buffers"
	cd protobuf-3.4.1
  mkdir -p $TOOLS_DIR/protobuf
	./configure --prefix=$TOOLS_DIR/protobuf CC=$C_COMPILER CXX=$CXX_COMPILER
	make -j 8
	make install
	echo "Completed"
	echo
	cd ../../
fi
export PATH=$TOOLS_DIR/protobuf/bin:$PATH
export LD_LIBRARY_PATH=$TOOLS_DIR/protobuf/lib:$LD_LIBRARY_PATH

# Install or Load Boost (with MPI and serialization).
if [ -f "$TOOLS_DIR/boost/lib/libboost_mpi.so" ] || [ -f "$TOOLS_DIR/boost/lib/libboost_mpi.dylib" ]; then
	echo "Found cached Boost"
else
	echo "Downloading Boost"
  mkdir -p downloads
  cd downloads
	wget -O boost_1_65_1.tar.bz2 https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2
	tar xjf boost_1_65_1.tar.bz2
	echo "Configuring and building Boost"
	cd boost_1_65_1
  mkdir -p $TOOLS_DIR/boost
  ./bootstrap.sh
  echo 'libraries =  --with-mpi --with-serialization ;' >> project-config.jam
  echo 'using mpi : mpic++ ;' >> project-config.jam
	echo 'using gcc : 6 ;' >> project-config.jam
	./b2 -j8 --prefix=$TOOLS_DIR/boost install
	echo "Completed"
	echo
	cd ../../
fi
export PATH=$TOOLS_DIR/boost/bin:$PATH
export LD_LIBRARY_PATH=$TOOLS_DIR/boost/lib:$LD_LIBRARY_PATH

# Install or Load Lockless.
if [ -f "$TOOLS_DIR/lockless/lib/libllalloc.so" ]; then
	echo "Found cached Lockless"
else
	echo "Downloading Lockless"
  mkdir -p downloads
  cd downloads
	wget -O lockless-1.3.tar.gz https://github.com/jl2922/lockless_allocator/archive/lockless-1.3.tar.gz
	tar xzf lockless-1.3.tar.gz
	echo "Configuring and building Lockless"
	cd lockless_allocator-lockless-1.3
  mkdir -p $TOOLS_DIR/lockless/lib
	make -j 8
	cp libllalloc.so.1.3 $TOOLS_DIR/lockless/lib/
	ln -sf $TOOLS_DIR/lockless/lib/libllalloc.so.1.3 $TOOLS_DIR/lockless/lib/libllalloc.so.1
	ln -sf $TOOLS_DIR/lockless/lib/libllalloc.so.1.3 $TOOLS_DIR/lockless/lib/libllalloc.so
	echo "Completed"
	echo
	cd ../../
fi
export LD_LIBRARY_PATH=$TOOLS_DIR/lockless/lib:$LD_LIBRARY_PATH

# Download Eigen.
echo "Downloading Eigen"
wget -O eigen-3.3.4.tar.bz2 http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2
tar xjf eigen-3.3.4.tar.bz2
rm eigen-3.3.4.tar.bz2
mkdir -p $TOOLS_DIR/eigen/include
mv eigen-eigen-5a0156e40feb/Eigen $TOOLS_DIR/eigen/include/
echo "Completed"
echo

# Download Google Test.
echo "Downloading Google Test"
wget -O release-1.8.0.tar.gz https://github.com/google/googletest/archive/release-1.8.0.tar.gz
tar xzf release-1.8.0.tar.gz
rm release-1.8.0.tar.gz
mv googletest-release-1.8.0 gtest
echo "Completed"
echo

cp ci.mk local.mk
make proto
make -j
make all_tests -j