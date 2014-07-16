#Id: Makefile,v 1.2 2014-01-16 17:46:07-08 - - $

JAVASRC    = genbankFeature.java mtAnnotate.java
SOURCES    = ${JAVASRC} Makefile
ALLSOURCES = ${SOURCES}
MAINCLASS  = mtAnnotate
CLASSES    = ${patsubst %.java, %.class, ${JAVASRC}}
JARCLASSES = ${CLASSE}
JARFILE    = mtAnnotate
LISTING    = Listing.p
JARS	   = biojava3-alignment-3.0.8.jar biojava3-core-3.0.8.jar

all : ${JARFILE}

${JARFILE} : ${CLASSES}
	echo Main-class: ${MAINCLASS} >Manifest
	echo Class-Path: . biojava3-alignment-3.0.8.jar biojava3-core-3.0.8.jar >> Manifest
	jar cvfm ${JARFILE} Manifest ${JARCLASSES} ${JARS}
	chmod +x ${JARFILE}
	- rm Manifest

genbankFeature.class : genbankFeature.java
	javac $<

mtAnnotate.class : mtAnnotate.java
	javac -cp ./:./biojava3-alignment-3.0.8.jar:./biojava3-core-3.0.8.jar $<
clean :
	- rm ${JARCLASSES} Manifest

spotless : clean
	- rm ${JARFILE}
