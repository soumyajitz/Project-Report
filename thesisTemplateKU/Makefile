RM			= rm -f
PDFLATEX	= pdflatex
BIBTEX		= bibtex

all: thesis tidy

tidy:
	$(RM) *.log
	$(RM) *.aux
	$(RM) *.bbl
	$(RM) *.blg
	$(RM) *.dvi

clean: tidy
	$(RM) *.pdf
	$(RM) *.toc
	$(RM) *.lof
	$(RM) *.lot

show: clean all
	open thesis.pdf

thesis:	thesis.tex
	$(PDFLATEX) thesis
	$(BIBTEX) thesis
	$(PDFLATEX) thesis
	$(PDFLATEX) thesis
