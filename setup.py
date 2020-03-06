from setuptools import setup

setup(
    name='MCbuilder',
    version='0.1.0',
    description='A macrocomplex builder from pdb interaction files',
    url='https://github.com/anmeert/sbi-project.git',
    author='Annika, Elena and Paula',
    keywords='protein pdb macrocomplex',
    packages=[],
    scripts=[],
    license='',
    long_description=open('README.md').read(),
    install_requires=[
        "biopython >= 1.76 ",
        "numpy >= 1.18.1",
    ],
    python_requires='>=3',
)
