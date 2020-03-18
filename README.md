# sbi-project

## Index
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [Introduction](#Introduction)
- [Background and Scientific Explanation](#Background-and-Scientific-explanation)
- [Algorithm](#Algorithm)
- [Tutorial](#Tutorial)
- [Examples Analysis](#Examples-Analysis)
- [Limitations](#limitations)
- [References](#References)
<!-- /TOC -->

## Introduction

Due to the importance of knowing the structure of protein-protein interactions (PPIs) in the cell, the goal of this project is to perform an algorithm to modulate a protein macro-complex from individual pairs of interactions using Bioinformatics resources.

## Background and Scientific Explanation

Proteins are versatile molecules that play many critical roles in the body. Individual proteins are capable of producing a large variety of protein complexes, or even complexes with other molecules, such as DNA, that allow them to perform their function. Proteins with more than one polypeptide chain form complexes. These complexes are known as the quaternary structures of proteins. The importance of knowing the structure of these complexes is because it modulates the biological activity of the protein and the separation of the subunits often leads to the loss of functionality. It refers to the spatial arrangement of these chains as well as any interactions among them, be they covalent or non-covalent.

Experimentally determining the full structure of a protein complex or quaternary structure is very costly and time-consuming. Therefore, it is essential to develop computational techniques that combine the experimental information and the data obtained from high-throughput methods which can correctly model the protein complexes.

To be able to establish these macrocomplexes we require a starting point, composed of Protein DataBank (PDB) files which feature experimentally determined structures. In our case, these structures are the combination of chains within our complex. As mentioned, experimentally determining the structure of a complex is problematic and so a suitable alternative is determining the physical structure of sections of the complex, i.e. the interaction between two chains within these complex, be they homo or heterodimers.

The core of our program is based on superimposition, where if we have sufficient PDB files from our complex then we can superimpose these and build a much more complex quarternary structure, modelling any interactions which the experimentally obtained PDB files may not show.

Superimposition is based on aligning two protein chain such that their carbon backbone lie one on top of the other as much as possible, where one chain is fixed and the other is rotated and translocated to minimize their separation. Though secondary structure can easily be different between distant homologues, their tertiary structure is more conserved and so homology can be more easily ascertained from this level.

In this program, we establish superimposition to implement new proteins onto the macrocomplex. If two chains have a low RMSD (root-mean-square deviation) we can determine that they are the same chain and so interactions in the individual PDB files can be combined to give a greater image of the macrocomplex.

Moreover, an accurate superimposition does not imply that a third chain is in that position. To be able to confirm whether a certain chain is a part of a macrocomplex in that location it is also important to check its surrounding chains and whether there are any steric clashes. This is the case of chains where, even though the chains superimpose correctly, the chain with which itis interacting occupies the space of a protein that is already there. This means that the interacting chains should not be placed there in the complex.

[INCLUDE DIAGRAM FROM BIORENDER]

Finally, PDB files often only have a part of the prtotein strcture and so we can rebuild the complete quaternary structure by using the relevant amino acid sequences, from FASTA files, the build out our model of the structure.
## Algorithm
como todavia no tenemos bien diseñado el alg, no se bien que poner.
### Input

### Alignment
The first step, after obtaining the information contained in the input, was to perform a **pairwise alignment** for each chain using the module `pairwise2` from the package `Bio` of `Biopython`. The objective of this alignment is to obtain sequences with **more than 95% identity** in order to later perform a structural alignment between them.
