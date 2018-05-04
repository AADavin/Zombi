﻿
## Zombi 
##### A simulator of species, genes and genomes that accounts for extinct lineages.

----------

### **Introduction** ###

Zombi is a flexible platform of genome evolution that accounts for extinct lineages. This feature makes it especially interesting for those studying organisms in which there exists Lateral Gene Transfers, since transfer events normally take place between lineages that have left no surviving descendants.
However, it is possible to use it simply as a tool to create different a wide range of evolutionary scenarios. Zombi output was thought to be especially simple and easy to read, understand and parse. 
 
Zombi uses a Birth-death model to generate a species tree and then it simulates the evolution of genomes along this species tree. Genomes evolve undergoing events of duplication, transfer, loss, translocation and inversion. Each event can
affect a region of variable length in the genome. Finally, it is possible to simulate also the evolution of the sequences along the branches of the different gene trees.

Zombi can be of great interest to those who want to test different evolutionary hypothesis under simulations and need to use a fast and easy to use tool to generate species trees, gene trees or sequences.

Writen by Adrián A. Davín 

----------

Using ideas, suggestions and comments coming from: Théo Tricou, Thibault Latrille, Nicolas Lartillot, Vincent Daubin, Damiene de Vienne, Eric Tannier and Gergely J. Szollosi

Please, if you have any doubts or need a hand, **just contact me here: aaredav@gmail.com** I am very happy to help.
Same thing if you find any bug or something that should be fixed!

CC BY-SA 3.0

----------

### **Citation** ###

Please, if you use Zombi, cite:

 ** Still to fill **

### **Installation** ###

First, clone the repository to your computer

    git clone https://github.com/AADavin/ZOMBI

Second, you need **python 3.6** installed with **ETE3**, **numpy** and pyvolve

    pip3 install ete3 numpy pyvolve
        

### **Usage** ###

**There is a detailed Wiki that explains how Zombi works [here](https://github.com/AADavin/ZOMBI/wiki)** It takes around 15 minutes of your time reading it! But if you want to launch it
right away, just read this.

There are **three main modes** to run Zombi: **T** (species Tree), **G** (Genomes) and  **S** (Sequences) 

You must run the computations in sequential order. This means that:
Computing **genomes** requires having computed previously a species tree computed with the mode T. 
Computing **sequences** requires having computed previously **genomes** with the mode G.

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