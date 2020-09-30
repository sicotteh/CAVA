<!-- vim-markdown-toc GFM -->

* [CAVA v1.2.3 README](#cava-v123-readme)
    * [1 INTRODUCTION](#1-introduction)
    * [2 PUBLICATION](#2-publication)
    * [3 DEPENDENCIES](#3-dependencies)
    * [4 INSTALLATION ON LINUX OR MAC](#4-installation-on-linux-or-mac)
    * [5 RUNNING CAVA](#5-running-cava)
    * [6 LICENCE](#7-licence)

<!-- vim-markdown-toc -->

CAVA README
==================

1 INTRODUCTION
--------------

CAVA (Clinical Annotation of VAriants) is a lightweight, fast, flexible 
and easy-to-use Next Generation Sequencing (NGS) variant annotation tool. 
It implements a clinical sequencing nomenclature (CSN), a fixed variant 
annotation consistent with the principles of the Human Genome Variation 
Society (HGVS) guidelines, optimised for automated clinical variant 
annotation of NGS data. 

CAVA has been extensively tested on exome data and is being used in the 
Mainstreaming Cancer Genetics (MCG) programme which applies NGS to 
increase the availability and affordability of clinical testing of 
cancer predisposition genes.

*UPDATE* CAVA is no longer being supported by the Rahman Team. This fork
is used by the Mayo Clinic to annotate our research and clinical samples.


2 PUBLICATION
-------------

If you use CAVA, please cite:

Márton Münz, Elise Ruark, Anthony Renwick, Emma Ramsay, Matthew Clarke, 
Shazia Mahamdallie, Victoria Cloke, Sheila Seal, Ann Strydom, 
Gerton Lunter, Nazneen Rahman. CSN and CAVA: variant annotation tools 
for rapid, robust next-generation sequencing analysis in the clinical 
setting. Genome Medicine 7:76, doi:10.1186/s13073-015-0195-6 (2015).


3 DEPENDENCIES
--------------

To install and run CAVA you will need the following dependencies installed:
- Python 3
- GCC and GNU make
- virtualenv

If your system is missing GCC and GNU make, these can be installed as 

On a Mac, the easiest way to set them up is to install Xcode Command Line Tools.
On a Debian or Ubuntu Linux, they can be set up by installing the build-essential package:

`bash sudo apt-get install build-essential`

If not already installed on your system, virtualenv can be set up by:

`bash 
pip install virtualenv
`

4 INSTALLATION ON LINUX OR MAC
------------------------------

```bash 
git clone https://github.com/Steven-N-Hart/dicom_wsi.git
# optional to checkout release
# e.g. git checkout v.1.2.4
./install.sh
```

CAVA uses virtualenv and pip to manage all its extra dependencies, 
which means that it will not clutter up your system by installing 
things globally. Everything it installs will go into a sub-directory 
in the CAVA directory. If you delete CAVA then everything it has 
installed will also be deleted.

Once the installation script has finished successfully, CAVA is ready 
for use.


5 RUNNING CAVA
--------------

CAVA can be run with the following simple command:

```bash
cava -c config.txt -i input.vcf -o output
```

It requires three command line arguments: 
the name of the configuration file (-c), the name of the input file (-i) 
and the prefix of the output file name (-o). 

6 LICENCE
---------

CAVA is released under MIT licence (see the LICENCE file).

