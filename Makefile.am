## Process this file with automake to produce Makefile.in
## Author: Erik.Bartos@gmail.com

SUBDIRS = src Data Docs Hyperons Nucleons Pions Pions Transition
AUTOMAKE_OPTIONS = subdir-objects

buenavistadocdir = ${prefix}/doc/buenavista
buenavistadoc_DATA = \
	README\
	COPYING\
	AUTHORS\
	ChangeLog\
	INSTALL\
	NEWS

EXTRA_DIST = $(buenavistadoc_DATA) 

# Remove doc directory on uninstall
uninstall-local:
	-rm -r $(buenavistadocdir)
