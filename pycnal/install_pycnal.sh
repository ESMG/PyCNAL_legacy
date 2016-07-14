#!/bin/sh

#DESTDIR=/usr/local
DESTDIR=$HOME/python
PYCNAL_PATH=$HOME/python/lib/python3.5/site-packages/pycnal
CURDIR=`pwd`

echo
echo "installing pycnal..."
echo
python setup.py build
python setup.py install --prefix=$DESTDIR

echo "installing external libraries..."
echo "installing gridgen..."
cd $CURDIR/external/nn
./configure --prefix=$DESTDIR
make install
cd $CURDIR/external/csa
./configure --prefix=$DESTDIR
make install
cd $CURDIR/external/gridutils
./configure CPPFLAGS=-I$DESTDIR/include LDFLAGS=-L$DESTDIR/lib CFLAGS=-I$DESTDIR/include --prefix=$DESTDIR
make install
cd $CURDIR/external/gridgen
export SHLIBS=-L$DESTDIR/lib
./configure CPPFLAGS=-I$DESTDIR/include LDFLAGS=-L$DESTDIR/lib CFLAGS=-I$DESTDIR/include --prefix=$DESTDIR
make
make lib
make shlib
make install
#
#PYCNAL_PATH=`python -c 'import pycnal ; print(pycnal.__path__[0])'`
#
# Now setting this above because this gave me:
# $ echo $PYCNAL_PATH
# scrip.so not found. Remapping function will not be available
# /u1/uaf/kshedstrom/python/lib/python3.5/site-packages/pycnal
#
# One could fix it...
cp libgridgen.so $PYCNAL_PATH
echo "installing scrip..."
cd $CURDIR/external/scrip/source
perl -pe "s#\/usr\/local#$DESTDIR#" makefile > makefile2
make -f makefile2
make -f makefile2 f2py
make -f makefile2 install
# Write it this way for Darwin...
cp -r scrip*.so* $PYCNAL_PATH
cd $CURDIR
echo
echo "Done installing pycnal..."
echo "pycnal make use of the so-called gridid file to store"
echo "grid information like the path to the grid file, the"
echo "number of vertical level, the vertical transformation"
echo "use, ..."
echo "Please set the environment variable PYCNAL_GRIDID_FILE"
echo "to point to your gridid file. A gridid file template"
echo "is available here:"
echo "$LOCALDIR/python/pycnal/pycnal/gridid.txt"
read -p "Press any key to continue or Ctrl+C to quit this install"
echo
