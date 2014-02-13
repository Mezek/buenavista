## Process this file with automake to produce Makefile.in


bin_PROGRAMS =  \
	Pions

Pions_SOURCES =  \
	Pions.cpp

AM_CXXFLAGS = -L/usr/lib64/root \
	`root-config \
	--libs` \
	-lMinuit2 \
	-I/usr/include/root 
