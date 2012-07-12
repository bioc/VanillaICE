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
	R --vanilla --slave < ${SUITE}

trelease:
	R_LIBS_USER=~/Library/R/2.15-bioc-release/library R-2.15 CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/2.15-bioc-release/library R-2.15 --vanilla --slave < ${SUITE}

tdevel:
	R_LIBS_USER=~/Library/R/2.15-bioc-devel/library R-2.15 CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/2.15-bioc-devel/library R-2.15 --vanilla --slave < ${SUITE}

idevel:
	R_LIBS_USER=~/Library/R/2.15-bioc-devel/library R-2.15 CMD INSTALL --no-lock ${PKG}

irelease:
	R_LIBS_USER=~/Library/R/2.15-bioc-release/library R-2.15 CMD INSTALL --no-lock ${PKG}


