
## SimuLyon (working name) 
##### A new tool to simulate genome evolution accounting for dead lineages

----------

### **Introduction** ###

SimuLyon is a simulator of species tree and genome evolution that accounts for extinct lineages. 
This feature makes it especially interesting for those studying organisms in which there exists Lateral Gene Transfers, since transfer events can take place between lineages that have left no surviving descendants.
SimuLyon uses a Birth-death model to generate a species tree. Different functions are included to generate species tree with non-constant rates of speciation and extinction. 
Then it simulates the evolution of genomes along this species tree. Genomes evolve undergoing events of duplication, transfer, loss, translocation and inversion. 

SimuLyon can be of great interest to those who want to test different evolutionary hypothesis under simulations and need to use a fast and easy to use tool to generate species trees, gene trees and sequences.

Please read the manual before using it. I know you are in a hurry but it only takes around 15 minutes of your time. If you are really in a hurry
then just follow the **Example 1**.

**The current version of simuLyon runs in discrete time. A newer version in continous time will come soon**

Writen by Adrián A. Davín 
Using ideas, suggestions and comments coming from: (still to fill) 

contact to aaredav@gmail.com

----------

Watch out! simuLyon is **unfinished**. The current version can undergo big changes

========================================================================

### **Usage** ###

You need **python 3.6** installed with **ETE3** and **numpy**

There are **four** modes to run simuLyon: **T** (species Tree), **G** (Genomes), **F** (gene Families) and  **S** (Sequences) 

Computing **genomes** or **gene families** requires having computed previously a species tree computed with the mode T
Computing **sequences** requires having computed previously either **genomes** or **gene families**

The parameters are read from a tsv file that it can be modified with any text editor. 

To start using simuLyon right away you can use:

    python simuLyon.py T SpeciesTreeParameters.tsv /Output_folder

The Species Tree will be created in the /Output folder, along some other useful files (this is explained in the Output section of this manual)

Then, you can simulate the evolution of genomes in that species tree using:

    python simuLyon.py G GenomeParameters.tsv /previous_Output_folder

Make sure that the Output_folder is the same that you were using when you generated the Species Tree! 
simuLyon will simulate the evolution of genomes inside the whole species tree and it will print a very detailed
information about them

Go through the examples for more details    
   
### **Generating a species tree (T)** ###

SimuLyon uses a Birth-death model to generate a species tree. In each step of time, lineages can speciate given rise to two new lineages
or go extinct. Different functions are included to generate trees with variable rates, namely:

 - 0: Global rates of speciation and extinction - This is the "traditional" way to simulate species tree using a Birth and death model. Both rates are constant over time and over all lineages
 - 1: Lineage-specific rates (time and lineage autocorrelated)  - **Not ready yet**.
 - 2: Lineage-specific rates (lineage autocorrelated) - **Not ready yet**
 - 3: Lineage-specific rates (uncorrelated) - **Not ready yet**
 - 4: User defined rates - The user can modify  the speciation and extinction rates in different interval of times. See example 4 for more details

Once a species tree is generated some other files are generated along. Those files are used by simuLyon when computing genomes

### **Generating a genome (G)** ###

To generate a genome it is first necessary to simulate a Species Tree using simuLyon. And no, it is not possible to input an externally computed tree, **but it will be in future versions**

To simulate genomes, simuLyon starts with an ancestral genome at the root, with a given number of genes. For now in the current version, all genes families present in this 
genome have a single copy (so in this ancestral genome there are no duplicated genes). In the current version, genes have 4 fields, separated by underscored.
For example, a given gene can be:

n4_+_5_12

The meaning of this is:

* n4: The branch in the species tree where that gene is found
* +: The orientation of the gene (it can be + or -). This is determined randomly with a bernouilli distribution p=0.5
* 5: The gene family identifier (in this case, gene family 5)
* 12: An unique identifier for that gene. 
 
A genome is an ordered collection of genes. So if we begin with a genome that has 5 genes, what we see is something like

| Position | Gene_family | Orientation | Id |
|----------|-------------|-------------|----|
| 0        | 1           | +           | 1  |
| 1        | 2           | -           | 1  |
| 2        | 3           | +           | 1  |
| 3        | 4           | +           | 1  |
| 4        | 5           | -           | 1  |

The genome is circular, so the position 4 is adjacent to the position 3 and 0

The events of evolution in this mode are:

* D: Duplications. A segment of the genome is duplicated. The new copy is inserted next to the old one
* L: Losses. A segment of the genome is lost
* T: Transfers. A segment of the genome is transferred to a contemporary species. The segment is inserted in a random position
* C: Translocations. A segment of the genome changes its position within the genome
* I: Inversions. A segment of the genome inverts its position
* O: Originations. A new gene family appears and it is inserted in a random position

The rates in this case are **genome-wise**.

There is also an additional rate for each event. This is called the **extension_rate**. This number (between 0 and 1) is
the p parameters of a **geometric distribution** that controls the length of the affected segment.

For example, if DUPLICATION_E == 1, the extension of the segment duplicated will be always 1 (meaning that only one gene is duplicated at a time)

By changing this parameter we can fine tune the extension associated to the different events. If inversions affect normally large chunks of the genome,
it suffices to use a low p.

Origination of new gene families are always of size 1, meaning that it is not possible to have an origination of two gene families in the same step of time.

Once that the full evolution of genomes has been simulated, simuLyon prints also the gene trees associated to the different gene families, all the events taking place
in each gene family, the events taking place in each branch and the genomes of each node in the species tree.

### **Generating gene families (F)** ###

**Not working yet**

With this method is possible to simulate genes that evolve inside the species tree independent from other genes. Genes families evolve inside the species tree undergoing events 
of duplications, transfers and losses. This method however, does not provide any information regarding the position of the different genes within the genomes, but technically simulating
a large number of independent gene families it is possible to study the presence or absence of different genes in different positions of the species tree to recover genomes. In this case genomes are simply a bunch
of genes that are found in no particular order in a given point of the species tree. This method uses **gene-wise** rates.

### **Generating sequences (S)** ###

**This is not written yet**

### Examples 

#### Example 1 - Simulating a species tree and the genomes ####

First thing we do is we are going to change the parameters of the Species Tree.
 For that we create a new file using:
 
     cp SpeciesTreeParameters.tsv   example1_SpeciesTreeParameters.tsv
     
By default (STOPPING_RULE = 0), the tree evolves until time reaches TOTAL_TIME
Let us change this to obtain a tree with 6 living lineages. We open the file and we modify:
 
 * STOPPING_RULE = 3
 * N_LINEAGES = 6
 
To make the example reproducible, let us also change the seed.

* SEED = 237

Let us also use change the default and extinction rates. For that, we will use as the speciation rate: 3 * 10^4
and as extinction rate 1 * 10^-4

For that we make SPECIATION_P0 = 0.0003 and EXTINCTION_P0 = 0.0001. Once we have changed that, we can run the command
        
     python simuLyon.py T example1_SpeciesTreeParameters.tsv EXAMPLE_1
   
 Then we can go to the folder /EXAMPLE_1 and inspect the files that have been created there.
 
 The resulting tree has 6 extant lineages and it can be found in the file /EXAMPLE_1/WholeTree
 The surviving species can be found in /EXAMPLE1/ExtantTree 
 
Once that we have computed a species tree, we can simulate **the evolution of genomes** inside this species tree. We will change the default parameters.

     cp GenomeParameters.tsv  example1_GenomeParameters.tsv
        
First, we are going to use the SEED = 237 to have reproducible results.

Second, we will make STEM_FAMILIES = 100. This makes the genome of the root to be composed of 100 gene families ordered in a circle 

Third, we are going to make inversions and translocations to affect large sections of the genome. For that we will change:

TRANSLOCATION_E = 0.05
INVERSION_E = 0.05 

 Finally, we launch simuLyon using the command
  
     python simuLyon.py G example1_GenomeParameters.tsv EXAMPLE_1
      

### **Output** ###

##### Mode T

*WholeTree:* The whole species tree including the dead lineages, in newick format

*ExtantTree:* The surviving species tree, in newick format

*ParametersLog.tsv*: Parameters used to run the simulation

*SpeciesTreeEvents.tsv*: Events (speciation and extinction) taking place in the species tree

*LineagesInTime.tsv*: Lineages alive in each unit of time (**I have to change the format to make it more efficient**)

#### Mode G

*Complete genomes*: A folder with one file per node of the species tree. Each file contains information about the
genome composition.

*CompleteGenomes*: A folder with one file per node of the species tree. Each file contains information about the
genome composition.  

*EventsPerBranch*: A folder with one file per branch of the species tree. Each file contains information about the
events taking place in that branch. There are 4 fields. 

1. Time: The time at which the event takes place
2. Event: The type of event that takes place in a given time (Duplication D, Losses L, ArrivingTransfer AT, LeavingTransfer LT, Inversions I, Translocations C and Originations O)
3. Nodes: Some more information about the kind of event:

* D: Three numbers separated by semicolons, e.g. 14;15;16. Correspond to the identifier of the gene that duplicates, and the identifiers of the new two genes generated
* L: One number. Corresponds to the identifier of the gene that is lost
* AT: Five numbers. In order: The branch of the species tree where the transfer leaves, the identifier of the gene that is transferred, the identifier of the gene that remains in the donor genome, the branch of the recipient genome and the identifier of the gene that is inserted in the recipient genome
* LT: Exactly as AT
* I: One number. Identifier of the gene affected
* C: One number. Identifier of the gene affected
* O: One number. Identifier of the gene appearing (always 1)

4. Gene_family. Identifier of the gene family affected by the event

Please notice that in the case of events that affect to several genes, this will be reflected in the first column (several events taking place at the same unit of time)
 
 *EventsPerBranch*: A folder containing the gene trees corresponding to the evolution of the different families and the gene trees pruned so that only surviving genes are represented.
 Please notice that gene trees with fewer than 2 copies are simply not output
  
#### Mode F

**Not working yet**

*Transfers.tsv*: A complete list of the transfer events, with donor and recipients

#### Mode S 

**Not written yet**

### Parameters ###

You can change the parameters directly in the tsv files. Watch out, it is a tabular separated values file, so you do not want mess spaces and tabs. In the case you do not understand what a parameter does, you do not want to touch it.

#### Species tree evolution (T) ####

**SPECIES_EVOLUTION_MODE**

 - 0: Global rates of speciation and extinction (only working mode so far)

**SEED**

If the SEED is 0, the seed is randomly computed for each run. Any other number will be used as the input seed

**PRUNING**

If PRUNING == 1, the file ExtantTree containing only surviving lineages is generated

**STOPPING_RULE**

 - 0: Time stops arriving at TOTAL_TIME
 - 1: Tree evolves until a total of species = N_LINEAGES (extinct and alive) have been generated
 - 2: Tree evolves until a total of species = N_LINEAGES (extinct) have been generated
 - 3: Tree evolves until a total of species = N_LINEAGES (alive) have been generated

**TOTAL_TIME**

To use with the stopping rule 0

**N_LINEAGES**

To use in combination with the stopping rules 1,2 or 3

**SPECIATION_D, EXTINCTION_D**

Distribution probability for each type of event

 - **fixed**: In this case, the value of the corresponding parameter 0 will be use as the rate. For example, if you are planning to use a constant rate of 0.005 for speciations you must use: speciation_d = fixed ; speciation_p0 = 0.005
 - **uniform**: The value will be sampled from an uniform distribution between the corresponding values of p0 and p1. For example, if you want the speciation rate to be randomly sampled from an uniform distribution U(0.003,0.007) you should use speciation_d = uniform; speciation_p0 = 0.003; speciation_p1 = 0.007.
 - **normal**: The value will be sampled from a normal distribution between the corresponding values of p0 and p1. For example, if you want the speciation rate to be randomly sampled from a normal distribution N(0.003,0.007) you should use speciation_d = normal; speciation_p0 = 0.003; speciation_p1 = 0.007.

#### Genome evolution (G) ####

**MODE**

 - 0: Global rates (only working mode so far)
 
**PREFIX**

- Adds a prefix to the gene families

**DUPLICATION_D, TRANSFER_D, LOSS_D, INVERSION_D, TRANSLOCATION_D, ORIGINATION_D**

Distribution probability for each type of event (fixed, uniform or normal)

**DUPLICATION_E, TRANSFER_E, LOSS_E, INVERSION_E, TRANSLOCATION_E**

The value of the p parameter of a geometric distribution that determines the extension of the genome (measured in number of genes) 
affected for event

**REPLACEMENT_T**

Not doing anything yet

**STEM_FAMILIES**

Number of gene families present in the ancestral genome at the root

**MIN_GENOME_SIZE**

The minimal size for a given genome. Smaller genomes will not be affected by losses events


#### Gene family evolution (F) ####


**STOPPING_RULE**

 - 0: A fixed number of gene families will be simulated == N_FAMILIES. Notice that gene families introduced by STEM_FAMILIES are not accounted for in the above computation
   
**GENE_EVOLUTION_MODE**

 - 0: Global rates 

**DUPLICATION_D, TRANSFER_D, LOSS_D, INVERSION_D, TRANSLOCATION_D**

Distribution probability for each type of event

 - **fixed**: In this case, the value of the corresponding parameter 0 will be use as the rate. For example, if you are planning to use a constant rate of 0.5 for duplications you must use: duplication_d = fixed ; duplication_p0 = 0.5
 - **uniform**: The value will be sampled from an uniform distribution between the corresponding values of p0 and p1. For example, if you want the duplication rate is randomly sampled from an uniform distribution U(3,5) you should use duplication_d = uniform; duplication_p0 = 3; duplication_p1 = 5.
 - **normal**: The value will be sampled from a normal distribution between the corresponding values of p0 and p1. For example, if you want the duplication rate is randomly sampled from an normal distribution N(3,5) you should use duplication_d = normal; duplication_p0 = 3; duplication_p1 = 5.

**STEM_FAMILIES**

Number of gene families already in the stem clade

#### Sequence evolution ####

This is not written yet

### FROM NOW ON, OLD TEXT AND UNFINISHED PIECES. PLEASE JUST IGNORE

**SEQUENCE_EVOLUTION_MODE**

 - 0: No heterogeneity
 - 1: Lineage-specific heterogeneity (time and lineage autocorrelated)
 - 2: Lineage-specific heterogeneity (lineage autocorrelated)
 - 3: Lineage-specific heterogeneity (uncorrelated)
 - 4: User defined heterogeneity

Very frequently you will be interested in having a Species Tree that contains only extant species and potentially, only a fraction of all the surviving lineages.

For doing that, we resort to the script **sampLyon.py**. This script is going to prune and sample the different trees generated. For running it you can use:

    python sampLyon.py T /previous_Output_folder Sample_name N
    
T indicates the mode (prune and sample species Tree), then we input the folder we created before, we give a name to this sample and finally we introduce the fraction (N)
of species that will be sampled.

Once you have done that, if you are interested in using the gene trees you *have to* run sampLyon this way: 

    python sampLyon.py G /previous_Output_folder /previous_Sample_name

This will cut the gene trees removing extinct species and other species that were not including in the sample

#### Example 4 - Simulating a dataset with a massive extinction event NOT FINISHED YET####

In this example we will simulate a dataset in which a massive extinction event takes between
time 0.6 and time 0.7.

 First thing we do is we are going to change the parameters of the Species Tree.
 For that we create a new file using:
 
     cp SpeciesTreeParameters.tsv   example2_SpeciesTreeParameters.tsv
     
 Then we are going to open example1_SpeciesTreeParameters.tsv in any text editor and modify
 SPECIES_EVOLUTION_MODE to 4 (user defined rates).
 
 We are going to change also the speciation rate to a fixed number of 10 and the extinction rate to 3.
 For that we make SPECIATION_P0 = 10 and EXTINCTION_P0 = 3  
 
 Then we are going to use the file massive_extinction.tsv included in the folder SimuLYON/Example/1
 
 If we open it we see a single line, that includes in this order:
 
 The beginning of the rate change, the end of the rate change, the number by which the speciation rate is multiplied
 during that period and the number by which the extinction rate is multiplied during that period.
 
 In plain words, the file says that during the time 0.6 and until time 0.7, the speciation rate is maintained the same (multiplied by 1)
 and the extinction rate is multiplied by 15.
 
 We need to modify then the example1_SpeciesTreeParameters.tsv to include also USER_DEFINE_RATES = path_to_massive_extinction.tsv
 
 Once we have done that, we can run simuLyon by using:

     python simuLyon.py T example1_SpeciesTreeParameters.tsv EXAMPLE_2 

### To do

Species tree

 - 1: Lineage-specific rates (time and lineage autocorrelated)  - Not ready yet
 - 2: Lineage-specific rates (lineage autocorrelated) - Not ready yet
 - 3: Lineage-specific rates (uncorrelated) - Not ready yet
 - 4: User defined rates
 
 Gene families
 
  - 1: Family-wise rates
 - 2: Lineage-specific rates (time and lineage autocorrelated)
 - 3: Lineage-specific rates (lineage autocorrelated)
 - 4: Lineage-specific rates (uncorrelated)
 - 5: User defined rates
 
 
 **USER_DEFINED_RATES**

File containing the multipliers on the speciation and the extinction rate per unit of time, to be used when the species tree is simulated in mode 4
The format is (separated by tabs)
time_start	time_end	spec_mult	ext_mult