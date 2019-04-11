# hubViz
<!--
hubViz: A Novel Tool for Hub-centric Visualization
-->

hubViz (a novel tool for **hub**-centric **Vis**ualization) is a novel tool for hub-centric visualization which is based on a latent space joint model (LSJM) for visualization.
'hubviz' package provides computationally efficient and user friendly interface to fit the hubViz models. 
[The 'hubviz' vignette](https://github.com/dongjunchung/GPA/blob/master/inst/doc/GPA-example.pdf?raw=true) provides a good start point for the step-by-step data analysis using 'hubviz' package.The following help pages provide a good start point for the genetic analysis using the 'hubviz' package, including the overview of 'hubviz' package and the example command lines:

```
library(hubviz)
package?hubviz
class?hubviz
```

Installation
============ 

The stable versions of 'hubviz' package can be obtained from the following URLs:

Package source: [https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.tar.gz?raw=true](https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.tar.gz?raw=true)

Windows binary: [https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.zip?raw=true](https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.zip?raw=true)

Mac OS/X binary: [https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.tgz?raw=true](https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.tgz?raw=true)

To install the developmental versions of 'hubviz' package, it's easiest to use the 'devtools' package. Note that the ‘hubviz’ package depends on the ‘Rcpp’ package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("dongjunchung/hubviz")
```

References
==========
Jin IH, and Jeon M (2019), "GPA: A doubly latent space joint model for local item
and person dependence in the analysis of item response data," *Psychometrika*, 84: 236-260.
