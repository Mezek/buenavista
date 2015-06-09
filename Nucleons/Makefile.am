
bin_PROGRAMS = \
	Kompare\
	Nucleons\
	NucleonsD\
	SingleGraph\
	NucleonsM\
	NucleonsMC\
	MultiGraph

Nucleons_SOURCES = \
	Nucleons.cpp

NucleonsD_SOURCES =  \
	Nucleons-D.cpp

SingleGraph_SOURCES =  \
	SingleGraph.cpp

NucleonsM_SOURCES =  \
	NucleonsCS.cpp

NucleonsMC_SOURCES =  \
	NucleonsMC.cpp

MultiGraph_SOURCES =  \
	MultiGraph.cpp

Kompare_SOURCES = \
	Komparision.cpp

AM_CXXFLAGS = \
	-L/usr/lib64/root \
	`root-config \
	--libs` \
	-lMinuit2 \
	-I/usr/include/root

## File created by the gnome-build tools


SUBDIRS = \
	Iachello \
	Delta
