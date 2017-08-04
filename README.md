InGRiD: Integrative Genomics Robust iDentification of cancer subgroups
===

InGRiD is a statistical approach to improve prediction of cancer subgroups and identification of key genes and pathways by integrating information from biological pathway databases.

Installation
===========

To install the development version of the INGRID package, it's easiest to use the 'devtools' package.

```
#install.packages("devtools")
library(devtools)
install_github("dongjunchung/INGRID")
```

Usage
===========

The R package vignette will provide a good start point for the genetic analysis using the INGRID package, including the overview of INGRID package and the example command lines:

```
library(INGRID)
vignette("INGRID-example")
```
The following two help pages will also provide quick references for INGRID package and the example command lines:

```
package?INGRID
class?INGRID
```

References
==========

Wei W, Sun Z, da Silveira WA, Yu Z, Lawson A, Hardiman G, Kelemen LE, and Chung D (2017), "Semi-supervised Identification of Cancer Subgroups using Survival Outcomes and Overlapping Grouping Information."
