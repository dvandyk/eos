CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@

lib_LTLIBRARIES = libeostaudecays.la
libeostaudecays_la_SOURCES = \
    tau-to-3-leptons.cc tau-to-3-leptons.hh \
	observables.cc observables.hh
libeostaudecays_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeostaudecays_la_LDFLAGS = $(GSL_LDFLAGS)
libeostaudecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	-lgslcblas \
	-lgsl

include_eos_taudecaysdir = $(includedir)/eos/tau-decays
include_eos_taudecays_HEADERS = \
	tau-to-3-leptons.hh \
	observables.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
    tau-to-3-leptons_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeostaudecays.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
tau_to_3_leptons_TEST_SOURCES = tau-to-3-leptons_TEST.cc