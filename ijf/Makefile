clean:
	rm -f *.dvi *.log *.aux *.spl *.bbl *.blg *.bak

build: clean
	pdflatex ForecastHybrid
	bibtex ForecastHybrid
	pdflatex ForecastHybrid
	pdflatex ForecastHybrid
