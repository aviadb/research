fname=main
${fname}.pdf: ${fname}.tex \
	General/abstract.tex \
	General/acknowledgments.tex \
	General/committee.tex \
	ANN/ann.tex \
	ANN/ann.bib \
	Features/features.tex
	if [ -e ${fname}.aux ]; \
	then \
	rm ${fname}.aux; \
	fi;
	pdflatex ${fname}
	bibtex ${fname}
	bibtex ${fname}1-blx
	bibtex ${fname}2-blx
	bibtex ${fname}3-blx
	# Add more if you have more chapters
	pdflatex ${fname}
clean:
	rm -f *.aux
	rm -f */*.aux
	rm -f *.bbl
	rm -f *.blg
	rm -f *.fls
	rm -f *.lof
	rm -f *.log
	rm -f *.lot
	rm -f *.toc
	rm -f *-blx.bib
	rm -f *.out
	rm -f *.run.xml
	rm -f pdfa.xmpi
	rm -f *.fdb_latexmk
	rm -f '${fname}.synctex(busy)'
	rm -f *.synctex.gz
	rm -f ${fname}.pdf