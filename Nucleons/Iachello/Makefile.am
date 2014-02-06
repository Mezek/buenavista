## Process this file with automake to produce Makefile.in


bin_PROGRAMS =  \
	NucIachello

NucIachello_SOURCES =  \
	NucIachello.cpp

AM_CXXFLAGS = -L/usr/lib64/root \
	`root-config \
	--libs` \
	-lMinuit2 \
	-I/usr/include/root 
