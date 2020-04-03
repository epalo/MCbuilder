import setuptools

setuptools.setup(
    name='MCbuilder',
    version='0.1.0',
    description='A macrocomplex builder from pdb interaction files',
    url='https://github.com/anmeert/MCbuilder',
    author='Annika, Elena and Paula',
    keywords='protein pdb macrocomplex',
    packages=['src'],
    scripts=['src/macrocomplex_builder.py','src/InteractingChain.py','src/Interaction.py','src/Complex.py','src/UserInteraction.py'], 
    license='',
    long_description=open('README.md').read(),
    install_requires=[
        "biopython >= 1.76 ",
        "numpy >= 1.18.1"
    ],
    python_requires='>=3',
)
