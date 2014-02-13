## Process this file with automake to produce Makefile.in
## Author: Erik.Bartos@gmail.com

SUBDIRS = src po Data Docs Hyperons Nucleons Pions Pions

buenavistadocdir = ${prefix}/doc/buenavista
buenavistadoc_DATA = \
	README\
	COPYING\
	AUTHORS\
	ChangeLog\
	INSTALL\
	NEWS


INTLTOOL_FILES = intltool-extract.in \
	intltool-merge.in \
	intltool-update.in

EXTRA_DIST = $(buenavistadoc_DATA) \
	$(INTLTOOL_FILES)

DISTCLEANFILES = intltool-extract \
	intltool-merge \
	intltool-update \
	po/.intltool-merge-cache


# Remove doc directory on uninstall
uninstall-local:
	-rm -r $(buenavistadocdir)
