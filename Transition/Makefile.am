## Process this file with automake to produce Makefile.in


AM_CXXFLAGS = -L/usr/lib64/root \
	`root-config \
	--libs` \
	-lMinuit2 \
	-I/usr/include/root 

bin_PROGRAMS =  \
	plotEta \
	fitEtaPrime \
	plotEtaPrime \
	plotPiZero \
	fitEta \
	fitPiZero

plotEta_SOURCES =  \
	plotEta.cpp

fitEtaPrime_SOURCES =  \
	fitEtaPrime.cpp

plotEtaPrime_SOURCES =  \
	plotEtaPrime.cpp

plotPiZero_SOURCES =  \
	plotPiZero.cpp

fitEta_SOURCES =  \
	fitEta.cpp

fitPiZero_SOURCES =  \
	fitPiZero.cpp
