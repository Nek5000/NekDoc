TEX_FILES  =  $(wildcard *.tex)

all: Nek_doc Nek_users

Nek_doc: Nek_doc.pdf

Nek_users: Nek_users.pdf

Nek_doc.pdf: $(TEX_FILES)
	pdflatex -shell-escape -draftmode Nek_doc.tex
	pdflatex -shell-escape Nek_doc.tex

Nek_users.pdf: $(TEX_FILES)
	pdflatex -shell-escape -draftmode Nek_users.tex
	pdflatex -shell-escape Nek_users.tex
	
clean:
	rm -f *~ *.ilg *bak *.idx *.ind *.aux *.toc *.ps *.log *.lof *.loa
	rm -f *.bbl *.blg *.dvi *.out Nek_doc.pdf Nek_users.pdf *.ps  *.los *.lot *.tdo

.PHONY: clean 
