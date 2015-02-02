#-----------------------------------------------------------------------------
# KMlocal: A testbed for k-means clustering algorithms based on local search
# Version: 1.7
#-----------------------------------------------------------------------------
# Copyright (c) 2004-2010 David M. Mount and the University of Maryland
# All Rights Reserved.
#
# PLEASE READ THE FILE "Copyright.txt" FOR COPYRIGHT INFORMATION AND
# DISCLAIMER.
#-----------------------------------------------------------------------------
# Makefile for the k-means test and evaluation program
#
# Main target (makes the kmeans test program.
#
#	make			make everything
#	make document		make documentation
#	make validate		perform validation tests
#	make clean		delete temporary files
#	make realclean		delete temps and executables
#
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Basic definitions
#		BIN_DIR		where to store executables
#		SRC_DIR		kmltest source directory
#		TEST_DIR	validation tests directory
#		DOC_DIR		documentation directory
#-----------------------------------------------------------------------------

BIN_DIR		= bin
SRC_DIR		= src
TEST_DIR	= test
DOC_DIR		= doc

default: 
	cd $(SRC_DIR); make all; cd ..

validate: $(KMTEST)
	cd $(TEST_DIR); make

document:
	cd $(DOC_DIR); make

all:
	cd $(SRC_DIR); make all; cd ..
	cd $(TEST_DIR); make
	cd $(DOC_DIR); make
	

#-----------------------------------------------------------------------------
# Cleaning
#-----------------------------------------------------------------------------

clean:
	cd $(SRC_DIR); make clean; cd ..
	cd $(TEST_DIR); make clean; cd ..
	cd $(DOC_DIR); make clean; cd ..

realclean: 
	cd $(SRC_DIR); make realclean; cd ..
	cd $(TEST_DIR); make realclean; cd ..
	cd $(DOC_DIR); make realclean; cd ..
