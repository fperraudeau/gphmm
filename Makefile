#
# Generic R Package Makefile
# 
SHELL := /bin/bash -e

.SECONDEXPANSION:

prefix   ?= /usr/local
bindir   ?= $(prefix)/bin
libdir   ?= $(prefix)/lib
etcdir   ?= $(prefix)/etc
sharedir ?= $(prefix)/share

ENV_VARS := prefix=$(prefix) bindir=$(bindir) libdir=$(libdir) \
	etcdir=$(etcdir) sharedir=$(sharedir)
MAKE_ENV := $(MAKE) $(ENV_VARS)

## R-package related.
CRAN := 'http://cran.us.r-project.org'
BIOC := 'http://bioconductor.org/biocLite.R'

## Packages used to install/test.
R_PREREQS = devtools testthat roxygen2 docopt 

R_LIBS ?= $(shell Rscript -e "cat(.libPaths()[1])")
RSCRIPT := export R_LIBS=$(R_LIBS) && Rscript --vanilla -e
PKG_NAME := $(notdir $(CURDIR))

env:
	@echo ----------------------------------------------------------------
	@echo Installing Package: $(PKG_NAME)
	@echo 
	@echo prefix    : $(prefix)
	@echo bin       : $(bindir)
	@echo lib       : $(libdir)
	@echo etc       : $(etcdir)
	@echo share     : $(sharedir)
	@echo R_LIBS    : $(R_LIBS)
	@echo
	@echo ----------------------------------------------------------------

install: pkg-install
	for f in $(shell ls inst/gr-*); do \
		X=`echo $$f | sed 's/inst\///g' | sed 's/.R//g'` \
			&& cp $$f $(bindir)/$$X	&& chmod a+x $(bindir)/$$X; \
	done

## Really, all packages should use testthat, but as a holdover,
## tiny-test-harness is also supported.
TEST_CODE = $(if $(shell ls tests/tth.R),\
	"source('tth.R'); runTests('test-.*.R$$');",\
	"require(testthat); test_check('$(PKG_NAME)');")
test: install
	@echo Running tests in $(PKG_NAME)
	@cd tests && $(RSCRIPT) $(TEST_CODE) && cd - > /dev/null

clean:
	@for P in ".o$$" ".so$$" ".Rhistory$$" "~$$"; do \
		find . -type f | grep "$${P}" | xargs rm -f; done
uninstall:
	R CMD REMOVE $(PKG_NAME)
	for f in $(shell ls inst/gr-*); do X=`echo $$f | sed 's/inst\///g' | sed 's/.R//g'` \
		&& rm -f $(bindir)/$$X; done

roxygen: r-prereqs
	$(RSCRIPT) "require(roxygen2); roxygenize();"

## Rebuild if any file in package has modified time after the
## package's DESCRIPTION file in the installed location.
##
ALL_FILES = $(shell find $(CURDIR) -type f)
PKG_DESC = $(addsuffix /DESCRIPTION, $(addprefix $(R_LIBS)/, $(PKG_NAME)))
pkg-install: r-prereqs $(PKG_DESC)
$(PKG_DESC): $(ALL_FILES)
	$(RSCRIPT) "source($(BIOC)); devtools::install(repos = biocinstallRepos());"

## Install prereqs if not installed.
##
R_PREREQS_INDEX = $(addsuffix /INDEX, $(addprefix $(R_LIBS)/, $(R_PREREQS)))
r-prereqs: $(R_PREREQS_INDEX)
$(R_PREREQS_INDEX):
	$(eval PKGNAME := $(subst /INDEX,,$(subst $(R_LIBS)/,,$@)))
	$(RSCRIPT) "if (! require($(PKGNAME))) install.packages('$(PKGNAME)', repos = $(CRAN))"

.PHONY: install uninstall clean env test roxygen

##
## I'm torn as to whether this should really go in here, but it's
## probably the only really sane way of managing the dependency.
pbsim-install: $(bindir)/pbsim
$(bindir)/pbsim: 
	mkdir -p pbsim-tmp
	cd pbsim-tmp && \
	wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/pbsim/pbsim-1.0.3.tar.gz && \
	tar zxvf pbsim-1.0.3.tar.gz && cd pbsim-1.0.3 && \
	./configure --prefix=$(prefix) && make && make install && \
	cd ../.. && rm -rf pbsim-temp
