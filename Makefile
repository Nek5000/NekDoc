TEX_FILES  =  $(wildcard *.tex)
BUILD_DIR  = build

all: pdf html

gh-pages: pdf html
	cd $(BUILD_DIR) && \
		git init && \
		git checkout -b gh-pages && \
		git add . && \
		git remote add origin git@github.com:Nek5000/NekDoc.git && \
		git commit -m "Updating Project Page" && \
		git push origin gh-pages --force

html: $(BUILD_DIR)/Nek_users.html

$(BUILD_DIR)/Nek_users.html: $(TEX_FILES) | $(BUILD_DIR)
	htlatex Nek_users.tex "html,index=1,2,fn-in" -cdefault -d$(BUILD_DIR)/

pdf: $(BUILD_DIR)/Nek_users.pdf

$(BUILD_DIR)/Nek_users.pdf: $(TEX_FILES) | $(BUILD_DIR)
	pdflatex -shell-escape -draftmode Nek_users.tex
	bibtex Nek_users
	pdflatex -shell-escape -draftmode Nek_users.tex
	pdflatex -shell-escape -output-directory $(BUILD_DIR) Nek_users.tex
	
clean:
	rm -f *~ *.ilg *bak *.idx *.ind *.aux *.toc *.ps *.log *.lof *.loa
	rm -f *.bbl *.blg *.dvi *.out Nek_users.pdf *.ps  *.los *.lot *.tdo
	rm -f *.html *.css *.4ct *.4tc *.idv *.lg *.tmp *.xref
	rm -rf $(BUILD_DIR)

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

.PHONY: clean 
