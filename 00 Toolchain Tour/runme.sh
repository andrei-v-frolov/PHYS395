#!/bin/bash

gfortran -O3 -fdefault-real-8 hello.f90 -o hello
./hello > DATA
