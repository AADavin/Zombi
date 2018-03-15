
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

You need **python 3.6** installed with **ETE3** and **numpy**

    pip3 install ete3 numpy

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

Go through the examples for more details    
   
### **Generating a species tree (T)** ###

SimuLyon uses a Birth-death model to generate a species tree. Lineages can speciate given rise to two new lineages
or go extinct over time. 

### **Generating a genome (G)** ###

To generate a genome it is first necessary to simulate a Species Tree using simuLyon. And no, it is not possible to input an externally computed tree, **but it will be in future versions**

To simulate genomes, simuLyon starts with an ancestral genome at the root, with a given number of genes. For now in the current version, all genes families present in this 
genome have a single copy (so in this ancestral genome there are no duplicated genes). 

A genome is an ordered collection of genes. So if we begin with a genome that has 5 genes, what we see is something like


| Position | Gene_family | Orientation | Id |
|----------|-------------|-------------|----|
| 0        | 1           | +           | 1  |
| 1        | 2           | -           | 1  |
| 2        | 3           | +           | 1  |
| 3        | 4           | +           | 1  |
| 4        | 5           | -           | 1  |


The meaning of this is:
* Position: The position in the genome. The genome is circular, so the position 4 is adjacent to the position 3 and 0
* Gene_family: The identifier of the gene family
* Orientation: The orientation of that gene in the genome
* Id: The identifier of the gene.

Genomes evolve undergoing a series of events:

* **D:** **Duplications**. A segment of the genome is duplicated. The new copy is inserted next to the old one
* **L:** **Losses**. A segment of the genome is lost
* **T:** **Transfers**. A segment of the genome is transferred to a contemporary species. The segment is inserted in a random position. Transfers can be replacement transfers 
* **C:** **Translocations**. A segment of the genome changes its position within the genome
* **I:** **Inversions**. A segment of the genome inverts its position
* **O:** **Originations**. A new gene family appears and it is inserted in a random position

The rates in this case are **genome-wise**. For instance, a duplication rate of 3 means 3 duplication events per genome per unit of time.

There is also an additional rate for each event. This is called the **extension_rate**. This number (between 0 and 1) is
the p parameters of a **geometric distribution** that controls the length of the affected segment.

For example, if DUPLICATION_EXTENSION == 1, the extension of the segment duplicated will be always 1 (meaning that only one gene is duplicated at a time)

By changing this parameter we can fine tune the extension associated to the different events. If inversions affect normally large chunks of the genome,
it suffices to use a low p.

Origination of new gene families are always of size 1, meaning that it is not possible to have an origination of two gene families in the same step of time.
Once that the full evolution of genomes has been simulated, simuLyon prints also the gene trees associated to the different gene families, all the events taking place
in each gene family, the events taking place in each branch and the genomes of each node in the species tree.

There are two other events that do not depend intrinsically on genomes but in the species tree that is used to simulate genome evolution

* **S**: **Speciation**. When a genome arrives at a speciation node, the genome is divided and continues to evolve in both descendant branches
* **E**: **Extinction**. When a genome arrives at a extinction event, the genome stop its evolution

----------

*Some advances details regarding the genes identifiers: You might want to skip this part if you are reading this for the first time*

Events that introduce nodes in the topology of the gene tree, change the identifier of the gene.
For example, let us say that in the root we have a gene whose identifier is 1. If the genome where the gene undergoes a speciation,
the two branches will inherit: one a gene whose identifier is 2 and the other one 3. A duplication will change also the identifiers of the
duplicated genes. When a gene has been transferred, it changes the identifier of the gene remaining in the genome and in the recipient genome.
This way is easy to track the events that have given rise to different tree topologies. Inversions and translocations do not introduce
changes in the tree topology and for that reason they do not change the identifier of the affected genes.

----------

### **Generating sequences (S)** ###

With this method it is possible to simulate an alignment for each gene that evolve inside the species tree.
This method can simulate non-coding DNA sequences, coding DNA sequences and amino-acids sequences.

This module requires Pyvolve (https://github.com/sjspielman/pyvolve), to install it use the command:

    pip3 install pyvolve
    
----------
    
## Advanced modes  

### Advanced modes for simulating species tree

#### Mode Tb - Branch-wise extinction and speciation rates

In this model, speciation and extinction rates are specific for each branch. Each time one new lineage emerges, a new value for
its extinction and speciation rates its sampled from a user defined distribution (normal or lognormal). The difference with the previous model
it is that rates are independently sample each time a new lineage emerge.

#### Mode Tp - Fine control of population

This model allows the user to fine control the size of the population

### Advanced modes for simulating genomes

#### Mode Gb - Branch-wise genome rates

In this model, genome event rates are specific for each branch of the species tree. Each time one new lineage emerges, a new value for
its genome events rates its sampled from a user defined distribution (normal or lognormal), in which the mean corresponds
to the value of the parent branch and the variance to the branch length

#### Mode Gs - Selection model genome rates

This model allows the user to take into account the importance of genes

#### Mode Gu - User control of genome rates

This model allows the user to fine control the genome rates

### Advanced modes for simulating sequences

#### Mode Sb - Branch-wise sequence rates
 
## Experimental modes  

The experimental modes are not ready to use yet

### Experimental modes for simulating species tree

#### Mode Ta - Autocorrelation of extinction and speciation rates

In this model, speciation and extinction rates are specific for each branch. Each time one new lineage emerges, a new value for
its extinction and speciation rates its sampled from a user defined distribution (normal or lognormal), in which the mean corresponds
to the value of the parent branch and the variance to the branch length

#### Mode Ti - Preparing SimuLyon with an input tree by the user 

This model allows the user to input a species tree

### Experimental modes for simulating genomes

#### Mode GT - Simultaneous evolution of species tree and genomes

This models simulates simulatenously the evolution of the species tree and the evolution of genomes

#### Mode Ga - Autocorrelation of genome rates

In this model, genome event rates are specific for each branch of the species tree. Each time one new lineage emerges, a new value for
its genome events rates its sampled from a user defined distribution (normal or lognormal), in which the mean corresponds
to the value of the parent branch and the variance to the branch length

### Advanced modes for simulating sequences

#### Mode Sa - Autocorrelated sequence rates
#### Mode Sf - Family-wise sequence rates
#### Mode Su - User control of sequence rates   


## **Output** ###

#### Mode T

**WholeTree.nwk** The whole species tree including the dead lineages, in newick format

**ExtantTree.nwk** The surviving species tree, in newick format

**Events.tsv**: Events (speciation and extinction) taking place in the species tree

#### Mode G

**Genomes**: A folder with one file per node of the species tree. Each file contains information about the
genome composition.  

**Gene_families**: A folder with one file per gene family. Each file contains information about the
events taking place in that gene family. There are 3 fields.  

 **1. Time**: The time at which the event takes place
 **2. Event**: The type of event that takes place in a given time (S, E, D, T, L, I, C, O and F. F stands for Final, meaning that the gene arrived *alive* till the end of the run)
 **3. Nodes**: Some more information about the kind of event:

 **S, D and T**: 6 fields separated by semicolons. This can be better understood looking at the picture:

[nodescode]: https://github.com/AADavin/SimuLyon/blob/master/Images/NodesCode.png "Nodes codes"
![alt text][nodescode]

* **L, I, C, O and F**: 2 fields separated by semicolons. First, the species tree branch where the event takes place and second, the identifier of the gene affected

**GeneTrees**: A folder containing the gene trees corresponding to the evolution of the different families and the gene trees pruned so that only surviving genes are represented.
 Please notice that gene trees with 2 or fewer surviving copies are simply not output. Mind that some families can appear at some point and then leave no surviving descendants!
 
 There are two types of trees:
 
 * **_wholetree.nwk**: A tree showing the whole evolution of that gene family
 * **_prunedtree.nwk**: A tree in which the genes that have not survived till the present time have been removed. **Normally you want to use this tree!**


**EventsPerBranch**: (Not output by default) A folder with one file per branch of the species tree. Each file contains information about the
events taking place in that branch. The codes are similar to the previously explained, but not the same. There are two main differences (for the sake of clarity). The first one is that transfers are divided into:

* **LT**: Leaving Transfers. Transfers that leave this branch
* **AT**: Arriving Transfers. Transfers that arrive at this branch.

The second difference is that the node of the nodes affected is given by:

GeneFamily_GeneIdentifier

So for example, if we go to the file n2_branchevents and we find the event L affecting at 4_3, means that the gene whose identifier is 3 belonging to the family 4 was lost in that branch in time given by the first column

Please also notice that in the case of events that affect to several genes, this will be reflected in the first column (several events taking place at the same unit of time)
  
   
 **Profiles**: (Not output by default) Here there is a file called Profiles.tsv that contains the node of the species tree in the columns and the gene families in the rows.
 The entries give the number of copies that each gene family has for each node of the species tree.
 
 #### Mode S

*Sequences*: A folder with one fasta file per gene of the species tree.
Each fasta alignment contains the simulated sequences obtained at the leaves of the tree, not the internal nodes.

### Examples (Not ready yet)
#### Example 1 - Simulating genomes for reconciliations
#### Example 2 - Simulating a massive extinction event
#### Example 3 - Simulating sequences for whatever

### Parameters ###

You can change the parameters directly in the tsv files. Watch out, it is a tabular separated values file, so you do not want mess spaces and tabs. In the case you do not understand what a parameter does, you do not want to touch it.

#### Species tree evolution (T) ####

**STOPPING_RULE**

 - 0: Time stops arriving at TOTAL_TIME
 - 1: Tree evolves until a total of species = N_LINEAGES (alive) have been generated
 
**TOTAL_TIME**

To use with the stopping rule 0

**N_LINEAGES**

To use in combination with the stopping rule 1. Otherwise is ignored

**SPECIATION, EXTINCTION**

Rates of evolution


#### Genome evolution (G) ####


**DUPLICATION, TRANSFER, LOSS, INVERSION, TRANSLOCATION, ORIGINATION**

The value for each type of event. 

**DUPLICATION_EXTENSION, TRANSFER_EXTENSION, LOSS_EXTENSION, INVERSION_EXTENSION, TRANSLOCATION_EXTENSION**

The value of the p parameter of a geometric distribution that determines the extension of the genome (measured in number of genes) 
affected for event

**REPLACEMENT_TRANSFER**

A number between 0 and 1 controling the probability of replacement transfers (they only happen if there is a homologous position in the recipient genome)

**STEM_FAMILIES**

Number of gene families present in the ancestral genome at the root

**MIN_GENOME_SIZE**

The minimal size for a given genome. Smaller genomes will not be affected by losses events

#### Sequence evolution (S) ####

