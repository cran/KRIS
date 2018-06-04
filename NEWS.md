---
title: "NEWS"
author: "Kridsadakorn Chaichoompu"
date: "04/06/2018"
output:
  html_document:
    keep_md: yes
    toc: yes
  pdf_document: 
    number_sections: yes
    toc: yes
---



# KRIS 1.1.1

## Updates

* Updated and added more examples in plot3view

# KRIS 1.1.0

## New features

* Added the function ```cal.pc.projection``` 

## Updates

* Moved ```replace.missing``` to for_parallel.R from principal.component.R

## Fixed bugs

* Fixed bugs in test_plot3views.R

# KRIS 1.0.4

## Updates

* Changed URL in the CITATION file from Gitlab to CRAN
* Changed email in the DESCRIPTION file from Gmail to Official email
* Fixed misspelling in a manual of data

---

# KRIS 1.0.3

## Updates

* Updated URLs and added details of package in the DESCRIPTION file
* Updated manual of all functions according from CRAN's comments

---

# KRIS 1.0.2

## Updates

* Fixed the help file for ```read.bed```, ```read.bed```, ```fst.each.snp.hudson``` and  ```fst.hudson```.

---

# KRIS 1.0.1

## Updates

* Fixed the help file for ```cal.pc.linear``` and  ```rubikclust```.

---

# KRIS 1.0.0.0

## Initial functions

* ```cal.pc.linear``` A function for linear principal component analysis (PCA)
* ```fst.each.snp.hudson``` A function for fixation index (Fst) calculation for 
all SNPs between two groups.
* ```fst.hudson``` A function for average fixation index (Fst) calculation 
between two groups.
* ```plot3views``` A function to create scatter plots in three views.
* ```read.bed``` Read the binary PLINK format (BED, BIM, and FAM)
* ```rubikclust``` A function for unsupervised clustering to detect rough 
structures and outliers.
* ```write.bed``` Write an list of SNP object to the binary PLINK format (BED, 
BIM, and FAM)
* ```xxt``` A function for calculating matrix multipication between a matrix and 
its transpose for large data.

## Initial R data 

* ```simsnp``` Synthetic dataset containing single nucleotide polymorphisms 
(SNP)
* ```sample_labels``` Synthetic dataset containing population labels for the 
dataset simsnp.

## Initial example files

* ```example_SNP.bed``` Synthetic dataset containing single nucleotide polymorphisms 
(SNP) in binary format
* ```example_SNP.bim``` Simulated SNP information
* ```example_SNP.fam``` Simulated sample information

## Updates

From the initial idea, some functions were changed their names:

* The name of function ```cal.PC.linear``` was changed to ```cal.pc.linear```.
* The name of function ```plot.3views``` was changed to ```plot3views```.
* The name of function ```rubikClust``` was changed to ```rubikclust```.
* The name of function ```XXt``` was changed to ```xxt```.
