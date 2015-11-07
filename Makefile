TEX_FILES  =  $(wildcard *.tex)

all: Nek_users

html: $(TEX_FILES)
	htlatex Nek_users.tex "html,index=1,2,fn-in"

Nek_users: Nek_users.pdf

Nek_users.pdf: $(TEX_FILES)
	pdflatex -shell-escape -draftmode Nek_users.tex
	pdflatex -shell-escape Nek_users.tex
	
clean:
	rm -f *~ *.ilg *bak *.idx *.ind *.aux *.toc *.ps *.log *.lof *.loa
	rm -f *.bbl *.blg *.dvi *.out Nek_users.pdf *.ps  *.los *.lot *.tdo

.PHONY: clean 
