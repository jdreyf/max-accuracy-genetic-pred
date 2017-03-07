# max-accuracy-genetic-pred
Code to accompany the paper: Dreyfuss JM, Levner D, Galagan JE, Church GM, Ramoni MF. How accurate can genetic predictions be? BMC Genomics. 2012 Jul 24;13:340. doi: 10.1186/1471-2164-13-340. PubMed PMID: 22827772; PubMed Central PMCID:
PMC3534619. http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-340.

To use these algorithms, call the functions: maxAUC, optSe, and/or optSp.

* maxAUC returns the maximal AUC given:
k=prevalence
pve=proportion of variance explained (e.g. heritability)
n.bins=number of bins (default: 100)

* optSe returns the optimal (maximal or minimal) sensitivity given:
k=prevalence
pve=proportion of variance explained (e.g. heritability)
sp=specificity
n.bins=number of bins (default: 1000)
direction="max" or "min" (default: "max")
thresh.vector=vector of thresholds in terms of bins (default: integer sequence from 1 to n.bins by 10)

* optSp returns the optimal (maximal or minimal) specificity given:
k=prevalence
pve=proportion of variance explained (e.g. heritability)
se=sensitivity
n.bins=number of bins (default: 1000)
direction="max" or "min" (default: "max")
thresh.vector=vector of thresholds in terms of bins (default: integer sequence from 1 to n.bins by 10)

* PARAMETERS  
k=prevalence. 0<=k<=1.
pve=proportion of variance explained (e.g. heritability). 0<=pve<=1.
se=sensitivity. 0<=se<=1.
sp=specificity.  0<=sp<=1.
n.bins=number of bins. Positive integer corresponding to ‘b’ in Detailed methods.
direction="max" or "min". Maximize or minimize sensitivity or specificity.
thresh.vector=vector of thresholds in terms of bins (default: integer sequence from 1 to n.bins by
10). This corresponds to ‘t’ in Detailed methods; people whose risk is >= threshold/n.bins test
positive. Thresholds should be positive integers <= n.bins.

* EXAMPLES  
Two examples follow. The commented-out numbers on the right-hand side are the numerical results I
got when using:
R version 2.12.1 (2010-12-16)
