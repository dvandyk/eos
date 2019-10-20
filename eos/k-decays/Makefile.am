CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@

lib_LTLIBRARIES = libeoskdecays.la
libeoskdecays_la_SOURCES = \
	k-to-l-nu.cc k-to-l-nu.hh
libeoskdecays_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeoskdecays_la_LDFLAGS = $(GSL_LDFLAGS)
libeoskdecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	-lgslcblas \
	-lgsl

include_eos_kdecaysdir = $(includedir)/eos/k-decays
include_eos_kdecays_HEADERS = \
	k-to-l-nu.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	k-to-l-nu_TEST

LDADD = \
	$(top_builddir)/test/libeostest.a \
	libeoskdecays.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
k_to_l_nu_TEST_SOURCES = k-to-l-nu_TEST.cc
