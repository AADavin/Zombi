

<img src="https://github.com/AADavin/Zombi/blob/master/Images/ZombiLogo.png" alt="zombilogo" height = "150" width="300"/>

### **Zombi: A phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages**

----------

### **Introduction** ###

**Zombi** is a flexible platform of genome evolution which can be of great interest to those who want to test different evolutionary hypotheses under simulations and need to use a fast and easy-to-use tool to generate species trees, gene trees or sequences.
Zombi's output is especially simple and easy to read, understand and parse. 


**Zombi current version**: **Zombi 1.0.0**


----------

Writen by Adrián A. Davín 

Using ideas, suggestions and comments coming from: Théo Tricou, Thibault Latrille, Nicolas Lartillot, Vincent Daubin, Damien de Vienne, Eric Tannier and Gergely J. Szollosi

Please, if you have any doubts or need a hand, **just contact me here: aaredav@gmail.com**


----------

### **Citation** ###

Please, if you use Zombi, cite:

**Zombi: A phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages.**

Adrian A. Davin, Theo Tricou, Eric Tannier, Damien M. de Vienne, Gergely J. Szollosi

Bioinformatics, btz710, https://doi.org/10.1093/bioinformatics/btz710

### **Installation** ###

First, clone the repository to your computer

    git clone https://github.com/AADavin/ZOMBI

Second, you need **python 3.6** installed with **ETE3**, **numpy**, **networkx** and **pyvolve**

    pip3 install ete3 numpy networkx pyvolve
    
You might find some problems when installing the ete3 library, more specifically, when ete3 tries to install the library six.
Simply use the command pip3 to install the missing libraries. For instance:

    pip3 install six
        

### **Usage** ###

**There is a detailed Wiki that explains how Zombi works [here](https://github.com/AADavin/ZOMBI/wiki)** and it takes around 15 minutes of your time reading it! But if you want to launch it
right away, just read this.  


There are **three main modes** to run Zombi: **T** (species Tree), **G** (Genomes) and  **S** (Sequences) 

You must run the computations in sequential order. This means that:
Computing **genomes** requires having computed previously a species tree computed with the T mode. 
Computing **sequences** requires having computed previously **genomes** with the G mode.


<img src="https://github.com/AADavin/Zombi/blob/master/Images/ZombiGraphicAbstract.png" alt="zombipipeline" height = "300" width="800"/>



The parameters are read from a .tsv file that can be modified with any text editor. 

To start using Zombi just write:

    python3 Zombi.py T ./Parameters/SpeciesTreeParameters.tsv ./Output_folder

These will generate a Species Tree along with some other files in ./Output folder/T

Then, you can simulate the evolution of genomes in that species tree using:

    python3 Zombi.py G ./Parameters/GenomeParameters.tsv ./Output_folder

Make sure that the Output_folder is the same that you were using when you generated the Species Tree! 
Zombi will simulate the evolution of genomes inside the whole species tree and it will print a very detailed
information about them in ./Output folder/G

Then, you can simulate the evolution of sequences for each gene in that species tree using:

    python3 Zombi.py S ./Parameters/SequenceParameters.tsv ./Output_folder
    
    
----------
    




**Zombi** is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

**Zombi** is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
<http://www.gnu.org/licenses/>.


