## Note that the tests are run on the scripts in the installed directly.
## Need to run 'make' to test any changes to the files in inst/unitTests
TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=doRUnit.R
#R=${R_HOME}/bin/R

all: inst test

inst: # Install package
	cd ${TOP}/..;\
	R CMD INSTALL ${PKG}

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.0-bioc-release/library R-3.0 --vanilla --slave < ${SUITE}

release: # Run unit tests
	cd ${TOP}/..;\
	R_LIBS_USER=~/Library/R/3.1-bioc-release/library R-3.1 CMD INSTALL ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.1-bioc-release/library R-3.1 --vanilla --slave < ${SUITE}

devel:
	R_LIBS_USER=~/Library/R/3.2-bioc-devel/library R-3.2 CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.2-bioc-devel/library R-3.2 --vanilla --slave < ${SUITE}
