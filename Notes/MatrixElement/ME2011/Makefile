# $Id: Makefile,v 1.1 2011/10/03 23:12:15 dmytro Exp $

DVIPS = dvips
DVIPSPREFLAGS = -t letter
DVIPSPOSTFLAGS = -Ppdf

# TEXINPUTS := $(HOME)/latex/foiltex:$(HOME)/tex-archive/macros/latex/contrib/supported/feynmf:.:$(TEXINPUTS)
# MPINPUTS := $(HOME)/tex-archive/macros/latex/contrib/supported/feynmf:.:$(MPINPUTS)

# LATEX = env TEXINPUTS=$(TEXINPUTS) latex
LATEX = pdflatex
LATEXFLAGS =

DVI = $(TEXSOURCES:.tex=.dvi)

SORCE = note.pdf 

ALLTEX = $(wildcard *.tex)
RERUN = '(No file .*\.toc|There were undefined references)'

.PHONY: all
all: $(SORCE)

%.pdf: %.tex $(wildcard *.tex) 
	@echo "====> LaTeX first pass"
	$(LATEX) $(LATEXFLAGS) $(<:.tex=)
	@if egrep -q $(RERUN) $*.log ; then echo "====> LaTeX rerun" && $(LATEX) $<; fi
	@if egrep -q $(RERUN) $*.log ; then echo "====> LaTeX rerun" && $(LATEX) $<; fi

.PHONY: clean
clean:
	rm -f *.pdf *.aux *.dvi *.log *.qsl *.sol *.lof *.lot *.toc *~
