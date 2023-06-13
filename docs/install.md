# System Requirements

COWSNPhR has been tested with Debian-based Linux systems, 
but should in principle work on any flavour of Linux, as well as MacOSX. 
Windows isn't supported, but it may very well be installable via bioconda. 

COWSNPhR should run on any regular desktop/laptop with 8 GB or RAM or more.

## Installing Using Conda

COWSNPhR is available as a conda package.


To install conda, as well as the COWSNPhR:

Download and run the conda installation scrips
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
conda update -q conda
```

Add the necessary channels to your conda installation
```
conda config --add channels olcbioinformatics
conda config --add channels conda-forge
conda config --add channels bioconda
```

Install the COWSNPhR conda package
```
conda install -c olcbioinformatics cowsnphr
```

With that done, typing `cowsnphr.py`, `cowsnphr` or `python -m cowsnphr_src.cowsnphr` launch the COWSNPhR pipeline. See the [Usage](usage.md) section for instructions on how to use the COWSNPhR.

## Run tests 

Test the installation by running the included tests

`git clone https://github.com/OLC-LOC-Bioinformatics/COWSNPhR.git`
`cd COWSNPhR`
`python -m pytest tests/test_vcf.py`
`python -m pytest tests/test_tree.py`