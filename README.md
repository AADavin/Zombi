
## SimuLyon (working name) 
##### A tool to simulate genome evolution accounting for dead lineages

----------

### **Introduction** ###

SimuLyon is a simulator of genome evolution that accounts for extinct lineages. 
This feature makes it especially interesting for those studying organisms in which there exists Lateral Gene Transfers, since transfer events can take place between lineages that have left no surviving descendants.
SimuLyon uses a Birth-death model to generate a species tree and then it simulates the evolution of genomes along this species tree. Genomes evolve undergoing events of duplication, transfer, loss, translocation and inversion. Each event can
affect a region of variable length in the genome. 

SimuLyon can be of great interest to those who want to test different evolutionary hypothesis under simulations and need to use a fast and easy to use tool to generate species trees, gene trees and sequences.

Please,**read the manual before using it!** 

Writen by Adrián A. Davín 

Using ideas, suggestions and comments coming from: (still to fill) 

**contact to aaredav@gmail.com**

Watch out! simuLyon is **unfinished**. The current version can undergo big changes

----------

### **Usage** ###

You need **python 3.6** installed with **ETE3**, **numpy** and pyvolve

    pip3 install ete3 numpy pyvolve

There are **three main modes** to run simuLyon: **T** (species Tree), **G** (Genomes) and  **S** (Sequences) 

You must run the computations in sequential order. This means that:
Computing **genomes** requires having computed previously a species tree computed with the mode T. 
Computing **sequences** requires having computed previously **genomes** with the mode G.

Each main mode has different advanced options. For example, it is possible to account for the variability in the rates of genome
evolution with the modes Ga or Gb. 

The parameters are read from a .tsv file that it can be modified with any text editor. 

To start using simuLyon right away you can use:

    python3 simuLyon.py T ./Parameters/SpeciesTreeParameters.tsv ./Output_folder

The Species Tree will be created in ./Output folder/T, along some other useful files (this is explained in the Output section of this manual)

Then, you can simulate the evolution of genomes in that species tree using:

    python3 simuLyon.py G ./Parameters/GenomeParameters.tsv ./Output_folder

Make sure that the Output_folder is the same that you were using when you generated the Species Tree! 
simuLyon will simulate the evolution of genomes inside the whole species tree and it will print a very detailed
information about them in ./Output folder/G

Then, you can simulate the evolution of sequences for each gene in that species tree using:

    python3 simuLyon.py S ./Parameters/SequenceParameters.tsv ./Output_folder