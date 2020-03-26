from setuptools import setup

setup(
    name='MCbuilder',
    version='0.1.0',
    description='A macrocomplex builder from pdb interaction files',
    url='https://github.com/anmeert/sbi-project',
    author='Annika, Elena and Paula',
    keywords='protein pdb macrocomplex',
    packages=['macrocomplex_builder'], # Creo que es como una carpeta donde tienes todos los scripts
    scripts=['macrocomplex_builder.py'], # Creo que aqui hay que poner los scripts necesarios para correr el macrocomplex_builder
    license='',
    long_description=open('README.md').read(),
    install_requires=[
        "biopython >= 1.76 ",
        "numpy >= 1.18.1"
    ],
    python_requires='>=3',
)
