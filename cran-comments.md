## Test environments
* local OS X install, R 3.4.1
* ubuntu 14.04 (on travis-ci), R 3.4.1
* win-builder (devel and release)

## R CMD check results
On local OS X, there were no ERRORs, WARNINGs, or NOTEs.

With win-builder, there were no ERRORs or WARNINGS, but one NOTE. The note is the following

checking CRAN incoming feasibility ... NOTE
Maintainer: 'Fanny Perraudeau <perraudeau.f@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  GPHMM (3:36, 4:29)
  
My email address is correct and I would like to keep GPHMM in capital letters.

## Downstream dependencies
There are currently no downstream dependencies for this package.
