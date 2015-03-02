CLEANFILES = \
	*~ \
	manual.acn \
	manual.acr \
	manual.alg \
	manual.aux \
	manual.bbl \
	manual.bcf \
	manual.blg \
	manual.glg \
	manual.glo \
	manual.gls \
	manual.ist \
	manual.out \
	manual.pdf \
	manual.run.xml \
	manual.tdo \
	manual.toc \
	manual.log
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = figures .

TEXSOURCES = \
	manual.tex \
	glossaries.tex \
	doc/installation.tex \
	doc/usage.tex \
	physics/b-baryon-hadronic-matrix-elements.tex \
	physics/b-meson-hadronic-matrix-elements.tex \
	physics/eft-b-to-s.tex \
	physics/incl-b-to-s-gamma.tex \
	physics/incl-b-to-s-dilepton.tex

BIBSOURCES = \
	references.bib

SOURCES = \
	$(BIBSOURCES) $(TEXSOURCES)

manual-am: manual.pdf

manual.pdf: $(SOURCES)
	pdflatex manual.tex
	pdflatex manual.tex
	makeglossaries manual
	biber manual
	pdflatex manual.tex
