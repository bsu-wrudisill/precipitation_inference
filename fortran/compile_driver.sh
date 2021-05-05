#!/bin/bash
gfortran -c snowmodule17.f90
gfortran -c driver.f90
f2py -c snowmodule17.f90 driver.f90 -m driver
