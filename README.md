# sbi-project

## Index
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [Introduction](#Introduction)
- [Background and Scientific explanation](#Background-and-Scientific-explanation)
- [Algorithm](#Algorithm)
- [Tutorial](#Tutorial)
- [Examples Analysis](#Examples-Analysis)
- [Limitations](#limitations) 
- [References](#References)
<!-- /TOC -->

## Introduction

Due to the importance of knowing the structure of protein-protein interactions (PPIs) in the cell, the goal of this proyect is to perform an algorithm to modulate a protein macro-complex from individual pairs of interactions using Bioinformatics resources. 

## Background and Scientific explanation

Proteins are versatil molecules that play many critical roles in the body. Individual proteins are capable of producing large variety of protein complexes, or even complexes with other molecules, such as DNA, that allow them to perform their function. Those complexes are known as quaternary structures of proteins. The importance of knowing the structure of these complexes is because of it modulates the biological activity of the protein and the separation of the subunits often leads to the loss of functionality. 

Experimentally determine the full structure of a protein complex or quaternary structure is very costly and time-consuming. Therefore, it is essential the development of computational techniques that combining the experimental information and the data obtained from high-througtput are able to correctly model the protein complexes. 



## Algorithm 
como todavia no tenemos bien diseñado el alg, no se bien que poner. 
### Input

### Alignment
The first step, after obtaining the information contained in the input, was to perform a **pairwise alignment** for each chain using the module `pairwise2` from the package `Bio` of `Biopython`. The objective of this alignment is to obtain sequences with **more than 95% identity** in order to later perform a structural alignment between them.