TXT  = $(shell ls *.txt | grep -v bcftools.txt)
HTML = $(TXT:.txt=.html)
XML  = $(TXT:.txt=.xml)

all: $(HTML)

%.html: %.txt header.inc footer.inc
	asciidoctor -b html5 -a linkcss -a stylesheet=index.css -a no-header-footer $<

clean:
	@rm -f $(HTML) $(XML)

