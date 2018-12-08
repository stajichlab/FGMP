#!/bin/bash

mkdir -p $PREFIX/src
mkdir -p $PREFIX/utils
mkdir -p $PREFIX/lib
mkdir -p $PREFIX/data
mkdir -p $PREFIX/sample

cp $SRC_DIR/src/* $PREFIX/src
cp -R $SRC_DIR/utils/* $PREFIX/utils
cp $SRC_DIR/README.md $PREFIX
cp $SRC_DIR/data/* $PREFIX/data
cp $SRC_DIR/fgmp.config $PREFIX
cp $SRC_DIR/lib/FGMP.pm $PREFIX/lib
cp $SRC_DIR/sample/* $PREFIX/sample
cp $SRC_DIR/sample_output/* $PREFIX/sample_output

chmod 755 $PREFIX/src/*
chmod 755 $PREFIX/utils/*
chmod 755 $PREFIX/lib/*
chmod 755 $PREFIX/bin

ln -s $PREFIX/src/fgmp.pl $PREFIX/fgmp
ln -s $PREFIX/lib/FGMP.pm $PREFIX/FGMP.pm
