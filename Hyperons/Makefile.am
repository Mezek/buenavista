## Process this file with automake to produce Makefile.in


bin_PROGRAMS =  \
	Mesons \
	Hyperons \
	NucleonsH

Mesons_SOURCES =  \
	Mesons.cpp

Hyperons_SOURCES =  \
	Hyperons.cpp

NucleonsH_SOURCES =  \
	NucleonsH.cpp

AM_CXXFLAGS = -L/usr/lib64/root \
	`root-config \
	--libs` \
	-lMinuit2 \
	-I/usr/include/root 
