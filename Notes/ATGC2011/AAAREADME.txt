======================================================================
                         How to get started
======================================================================
MacOS:
	tetex - has all that you need prepare the note
Linux:
	should work out of the box


======================================================================
                            Organization
======================================================================
LaTeX files:

        note.text - main file with little content in it. Please use a
	separate file for each major topic and include it in the main
	LaTeX file.

Plots:
	please use pdf format. ROOT by default may produce problematic
	pdf files, so storing output in eps and using epstopdf to get
	pictures may give you a better result

Style recommendations:
        
	http://cms-secr.web.cern.ch/cms-secr/Documents/cmspb/CMSGuidelines.pdf

Documentation for authors:

	http://cms-secr.web.cern.ch/cms-secr/Documents/cmspb/cmspb.html
	   
======================================================================
                            How To Run It
======================================================================

make

if you really want to do by hand, 

pdflatex note.tex
