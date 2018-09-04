

<img src="https://github.com/AADavin/Zombi/blob/master/Images/ZombiLogo.png" alt="zombilogo" height = "100" width="300"/>
### **A simulator of species, genes and genomes that accounts for extinct lineages**

----------

### **Introduction** ###

Zombi is a flexible platform of genome evolution which can be of great interest to those who want to test different evolutionary hypotheses under simulations and need to use a fast and easy to use tool to generate species trees, gene trees or sequences.
Zombi's output is especially simple and easy to read, understand and parse. 


----------

Writen by Adrián A. Davín 

Using ideas, suggestions and comments coming from: Théo Tricou, Thibault Latrille, Nicolas Lartillot, Vincent Daubin, Damien de Vienne, Eric Tannier and Gergely J. Szollosi

Please, if you have any doubts or need a hand, **just contact me here: aaredav@gmail.com**


----------

### **Citation** ###

Please, if you use Zombi, cite:

**Zombi: A simulator of species, genes and genomes that accounts for extinct lineages.**

Adrian A. Davin, Theo Tricou, Eric Tannier, Damien M. de Vienne, Gergely J. Szollosi

bioRxiv 339473; doi: https://doi.org/10.1101/339473

https://www.biorxiv.org/content/early/2018/06/07/339473

### **Installation** ###

First, clone the repository to your computer

    git clone https://github.com/AADavin/ZOMBI

Second, you need **python 3.6** installed with **ETE3**, **numpy** and pyvolve

    pip3 install ete3 numpy pyvolve
        

### **Usage** ###

**There is a detailed Wiki that explains how Zombi works [here](https://github.com/AADavin/ZOMBI/wiki)** and it takes around 15 minutes of your time reading it! But if you want to launch it
right away, just read this.

There are **three main modes** to run Zombi: **T** (species Tree), **G** (Genomes) and  **S** (Sequences) 

You must run the computations in sequential order. This means that:
Computing **genomes** requires having computed previously a species tree computed with the T mode. 
Computing **sequences** requires having computed previously **genomes** with the G mode.

The parameters are read from a .tsv file that can be modified with any text editor. 

To start using Zombi just write:

    python3 Zombi.py T ./Parameters/SpeciesTreeParameters.tsv ./Output_folder

These will generate a Species Tree along some other files in ./Output folder/T

Then, you can simulate the evolution of genomes in that species tree using:

    python3 Zombi.py G ./Parameters/GenomeParameters.tsv ./Output_folder

Make sure that the Output_folder is the same that you were using when you generated the Species Tree! 
Zombi will simulate the evolution of genomes inside the whole species tree and it will print a very detailed
information about them in ./Output folder/G

Then, you can simulate the evolution of sequences for each gene in that species tree using:

    python3 Zombi.py S ./Parameters/SequenceParameters.tsv ./Output_folder
    
    
----------
    

Distributed under:
CC BY-SA 3.0
