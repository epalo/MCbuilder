# Macrocomplex Builder
Macrocomplex Builder is a standalone `python3` application developed by **Annika Meert**, **Elena Pareja Lorente** and **Paula Torren Peraire**. The purpose is to model the macro-complex (quaternary) structure of biomolecules, formed by proteins and DNA or only proteins given a folder with the structures of protein-protein and protein-DNA interactions in PDB format. 
##Index


- [Requirements](#Requirements)
- [Installation](#Installation)
- [Package Structure](#Structure-of-the-Package)
- [Usage](#Usage)
- [Analysis of examples](#analysis-of-the-examples)
- [FAQS(?)](#FAQS)




## Requirements
For using Macrocomplex Builder, there are some requirements:

* [Python 3](https://www.python.org/)

**Third-party packages**

It was used several python modules for the tool development. These modules can be found in [PyPI](https://pypi.org/). 

* [biopython](https://pypi.org/project/biopython/)
* [numpy](https://pypi.org/project/numpy/)
<!--Preguntar si usamos algún paquete de estos más - no se si argparse y logging hay que ponerlo aquí. No me acuerdo si hacia falta instalarselo con PyPI-->


## Installation

To use the tool, the user first has to download the package from Git and then install it via `setup.py` in the python site-packages.

 ```bash
 $ git clone https://github.com/anmeert/sbi-project.git
 $ cd sbi-project
 $ sudo python3 setup.py install
 ``` 
## Structure of the Package
The package has been structured as follows:
<!--ASK ANNI, MAYBE THERE IS A BETTER WAY-->

![Tree](tree.png)

* `README.md`: 
* `report.md`:
* `UserInteraction.py`:
* `loggingSetup.py`:
* `macrocomplex_builder.py`:
* `processInputFiles.py`:
* `setup.py`:
* `get_example_files.py`:

## Usage

In this section, an explanation of how to run the script from the command line and what are the required and optional arguments is provided. 

The simples way to run it:

 ```bash
 $ macrocomplex_builder.py
 ```
In this case, it will run using the current folder as input. 

Use `-h` to get information about the rest of the arguments: 

```bash
 $ macrocomplex_builder.py -h 
# WHEN THE HELP IS DONE, COPY HERE 
# (-i, -o --> poner que sino pone nombre, se llamará macrocomplex.pdb modificar)
# (-chain_number, -stechiometry -v...)
 
```
## Analysis of the examples

## FAQS





