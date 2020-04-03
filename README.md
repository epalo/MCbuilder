# Macrocomplex Builder
Macrocomplex Builder is a standalone `python3` application developed by **Annika Meert**, **Elena Pareja Lorente** and **Paula Torren Peraire**. The purpose is to model the macro-complex (quaternary) structure of biomolecules, formed by proteins and nucleotids or only proteins given a folder with the structures of protein-protein and protein-DNA/RNA interactions in PDB format. 
##Index


- [Requirements](#Requirements)
- [Installation](#Installation)
- [Package Structure](#Structure-of-the-Package)
- [Tutorial](#Tutorial)
- [FAQS](#FAQS)




## Requirements
For using Macrocomplex Builder, there are some requirements:

* [Python 3](https://www.python.org/)

**Third-party packages**

It was used several python modules for the tool development. These modules can be found in [PyPI](https://pypi.org/). 

* [biopython](https://pypi.org/project/biopython/)
* [numpy](https://pypi.org/project/numpy/)

## Installation

To use the tool, the user first has to download the package from Git and then install it via `setup.py` in the python site-packages.

 ```bash
 $ git clone https://github.com/anmeert/sbi-project.git
 $ cd sbi-project
 $ sudo python3 setup.py install
 ``` 
## Structure of the Package
The package has been structured as follows:

![Tree](documentation/tree.png)

* **`Complex.py`**: A module used by `macrocomplex_builder.py` that contains the `Complex object` that stores a `Bio.PDB.Model.Model`, a list of `InteractingChains`, a `logging.getLogger(__name__) Object`, and (optionally) a dictionary for stoichiometry
* **`README.md`**: A MarkDown File containing tutorial and important information about the package. 
* **`report.md`**: A MarkDown File containing background and scientific explanation, analysis of the results and package limitations. 
* **`get_example_files.py`**: Complementary script for building the examples. 
* **`InteractingChain.py`**:A module used by `Complex.py` and `macrocomplex_builder.py` that contains the `InteractingChain Object` that stores a `Bio.PDB.Chain.Chain Object`, a file_index(int) that refers to the associated interaction-file, a string of residues as sequence and (optionally) the `InteractingChain` that its interacting with.
* **`Interaction.py`**: A module used by `macrocomplex_builder.py` that contains the `Interaction Object` stores a `Bio.PDB.Model.Model`, a first `InteractingChain` and a second `Interacting Chain`.
* **`UserInteraction.py`**: A module used by `Complex.py` and `macrocomplex_builder.py` that performs the import of user data and output of files.
* **`macrocomplex_builder.py`**: The main script to run the program. 
* **`setup.py`**: To install the package.


## Tutorial

In this section, an explanation of how to run the script from the command line and what are the required and optional arguments is provided. 

The simples way to run it:

 ```bash
 $ cd scr/
 $ macrocomplex_builder.py
 ```
In this case, it will run using the current folder as input and it will produce as an output a file named macrocomplex.pdb in the same folder. 

Input and output can also be specified (using `-i`and `-o` arguments). The input can be a folder with the interaction in pdb format and the sequences in fasta format or all the files to be used.
 
 ```bash
 $ macrocomplex_builder.py -i /path/pdb/fasta/folder -o /final/path/complex_name.pdb
 
 $ macrocomplex_builder.py -i /path/pdb/file.pdb /path/fasta/file.fasta -o /final/path/complex_name.pdb
 ```
Another option is to indicate the maximum number of strings that the user wants to add to the final complex. This number of strings may not be reached if it is not possible to add more strings due to stereo impediments. 

 ```bash
 $ macrocomplex_builder.py -i /path/pdb/fasta/folder -l 30
 ```
The user can also indicate the stoichiometry of the complex in a specific format:

 ```bash
 $ macrocomplex_builder.py -i /path/pdb/fasta/folder -s A3B5C8
 ```
It is also possible to change the default RSMD value before checking clashes. The default value is 0.5.

 ```bash
 $ macrocomplex_builder.py -i /path/pdb/fasta/folder -r 0.2
 ```

There are two run types: simple and complete. By default, it runs simply because it is less computationally expensive and gives quite good results. If the user would like to explore all the different possibilities the `-c` argument can be used (see report for more information about simple and complete algorithms). 

 ```bash
 $ macrocomplex_builder.py -i /path/pdb/fasta/folder -c
 ```
If the user wants to follow the program progression, verbose option `-v` will print the progression log and warnings in the standard output. By default, warnings are printed in the `log file`. 

 ```bash
 $ macrocomplex_builder.py -i /path/pdb/fasta/folder -v
 ```

Use `-h` to get information about all the arguments: 

```bash
 $ macrocomplex_builder.py -h 
 

usage: macrocomplex_builder.py [-h] [-i [INFILE [INFILE ...]]] [-o OUTFILE]
                               [-v] [-s STOICH] [-l LIMIT] [-c]
                               [-r RMSD_THRESHOLD]

Macrocomplex builder, creates protein macrocomplex from individual PDB files
and FASTA files.

optional arguments:
  -h, --help            show this help message and exit
  -i [INFILE [INFILE ...]], --input [INFILE [INFILE ...]]
                        Input fasta and PDB files, or directory containing
                        these
  -o OUTFILE, --output OUTFILE
                        PDB output file fo complete model
  -v, --verbose         Print progression log to standard output
  -s STOICH, --stoichiometry STOICH
                        Set stoichiometry for the macrocomplex. Input in
                        standard form (e.g. A1B4C6)
  -l LIMIT, --limit LIMIT
                        Limit number of chains in the protein
  -c, --complete        Run full recursion. Explores all possibilities
  -r RMSD_THRESHOLD, --rmsd RMSD_THRESHOLD
                        Set up rmsd threshold. Default is set to 0.5
```

## FAQS
* **What is the maximum number of chains that can be modelled?** The program can build models with no limit of chains, the problem is that it can only print 52 chains. Therefore, the maximum number of chains that the generated complex can contain is 52.
* **Can the algorithm handle with heteroatoms?** Yes, if the pdb inputs contain heteroatoms, these will be represented in the final complex.
* **Can the algorithm handle with different subunits?** Yes, as long as the coordinates of the interaction pairs are the same coordinates of the final complex. 










