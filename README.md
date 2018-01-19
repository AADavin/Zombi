
## SimuLyon  (working name)

----------

### **Introduction** ###

**SimuLyon** is a simulator of species tree and genome evolution that accounts for extinct lineages.

First, it uses a Birth-death model to generate a species tree.  Different functions are included to generate species tree with non constant rates of speciation and extinction

Second, genome evolution is simulated by the evolution of different gene families  along this species tree.  Gene families can be already present in the stem branch or appear with a given rate of origination along the different branches of the species tree. Genes can evolve through events of **duplication, loss and transfer**. Transfers can be or not **replacement transfers**.  SimuLyon also can account for the fact that **transfers occur preferentially between closely related lineages**. It also presents a detailed output of the genome evolution including the total number of families that have gone extinct and the ancestral genomes for the different inner nodes in the species tree.		

Third, the evolution of sequences (**TL**)

Please read the manual before using it. It should take around 15 minutes of your time.

Writen by Adrián A. Davín and ...

contact to aaredav@gmail.com

----------

========================================================================


### **Usage** ###

You need **python 3.6** installed with **ETE3** and **numpy**

Download the python scripts:

    simuLyon.py
    sampLyon.py

And the files:

    SpeciesTreeParameters.tsv
    GenomeParameters.tsv
    SequenceParameters.tsv

There are **three** modes to run simuLyon. **T** (species Tree) and **G** (Genomes) 
(and we will write a third one for sequences, **S**) -- **TL**

The order to obtain a complete dataset is 

 1. Species Tree
 2. Genomes
 3. Sequences

This means that if you want to get to the last point you must first compute the 2 first points. So, the very first thing you have to do is running simuLyon to create a species tree using simply

    python simuLyon.py T SpeciesTreeParameters.tsv /Output_folder

The Species Tree will be created in the /Output folder, along some other useful files (this is explained in the Output section of this manual)

Then, you can simulate the evolution of genomes in that species tree using:

    python simuLyon.py G GenomeParameters.tsv /Output_folder

Make sure that the Output_folder is the same!

These two steps generate a species tree (with extinct lineages). To generate the sequences, it is first recommended using sampLyon





Finally, we normally do not want to study the whole species tree but a smaller version of it. For that we use the script **sampLyon.py**. It cuts and prunes the Species Tree and the Gene Trees so that only a given fraction of the species surviving in present time are represented.

For example, if we generate a Species Tree with 10 000 species out of 1000 survive, but we are interested having about 10% of those surviving species in our final dataset, we can run the command

    python sampLyon.py 100 /Output_folder Sample100

sampLyon will cut the species tree to 100 randomly selected species (out of those living ones) and will output the whole information in the folder Sample100. Probably the best way to understand how SimuLyon and SampLyon work is by following the examples including in this manual.

### **Output** ###

####Mode S####


**WholeTree:** The whole species tree including the dead lineages, in newick format
**ParametersLog.tsv**: Parameters used to run the simulation
**SpeciesTreeEvents.tsv**: Events (speciation and extinction) taking place in the species tree
**LineagesInTime.tsv**: Lineages alive in each unit of time (I have to change the format to make it more efficient)


####Mode G ####

**Duplications.tsv**: Tsv including all the duplication taken place in each node of the whole tree
**Losses.tsv**: Tsv including all the losses taken place in each node of the whole tree
**LeavingTransfers.tsv**: Tsv including all the transfers leaving a given branch of the whole tree
**ArrivingTransfers.tsv**: Tsv including all the transfers arriving to a given branch of the whole tree
**Profiles.tsv**: Tsv including the number of copies present of each family for each node of the whole tree
**RawGeneFamilies**: Folder including all the gene trees of evolving inside the species tree

####Mode S ####

**TL**




### Examples ###

#### **Example 1 - Simulating a Massive Extinction event** ####

### Parameters ###

You can change the parameters directly in the tsv files. Watch out, it is a tabular separated values file, so you do not want mess spaces and tabs. In the case you do not understand what a parameter does, you do not want to touch it. 

#### Global parameters ####

This parameters are found in the python script and you are strongly encouraged not to change them

**TOTAL_TIME**
The total distance from the root of the species tree and the present time and the present. It is 1 by default. 

**TIME_INCREASE**

The size of the discretized time units. It is 0.0001 by default and recommended. A smaller number will produce more fine-grained scenarios at a larger computational cost. Think carefully before changing this.

#### Species evolution ####

**SPECIES_EVOLUTION_MODE**

 - 0: Global rates of speciation and extinction 
 - 1: Lineage-specific rates (time and lineage autocorrelated)
 - 2: Lineage-specific rates (lineage autocorrelated)
 - 3: Lineage-specific rates (uncorrelated)
 - 4: User defined rates

**SPECIATION_RATE**
 Speciation rate, measured in mean number of speciation per unit of time ( 

**EXTINCTION_RATE**
 Extinction rate, measured in mean number of speciation per unit of time   

**USER_DEFINED_RATES**
File containing the multipliers on the speciation and the extinction rate per unit of time, to be used when the species tree is simulated in mode 4
The format is (separated by tabs)
time_start	time_end	spec_mult	ext_mult	
See example 1 for more information

**MAX_SLEAVES **
The maximum number of leaves of the species trees at any instant. If you need to simulate very large trees you can modify this number. If you reach this number before the total age of the tree is TOTAL_TIME, the simulation is killed and another run starts. 


#### Genome evolution ####

**GENOME_EVOLUTION_MODE**

 - 0: A fixed number of gene families will be simulated == N_FAMILIES. Notice that gene families introduced by STEM_FAMILIES are not accounted for in the above computation
 - 1: Families will be simulated until the mean size of genomes >=  MEAN_SIZE_GENOME 
 **AA** NOT WRITTEN YET!!

**GENE_EVOLUTION_MODE**

 - 0: Global rates 
 - 1: Family-wise rates

Not written yet!

 - 2: Lineage-specific rates (time and lineage autocorrelated)
 - 3: Lineage-specific rates (lineage autocorrelated)
 - 4: Lineage-specific rates (uncorrelated)
 - 5: User defined rates

**DUPLICATION_D, TRANSFER_D, LOSS_D, REPLACEMENT_D**

Distribution probability for each type of event

 - **fixed**: In this case, the value of the corresponding parameter 0 will be use as the rate. For example, if you are planning to use a constant rate of 0.5 for duplications you must use: duplication_d = fixed ; duplication_p0 = 0.5
 - **uniform**: The value will be sampled from an uniform distribution between the corresponding values of p0 and p1. For example, if you want the duplication rate is randomly sampled from an uniform distribution U(3,5) you should use duplication_d = uniform; duplication_p0 = 3; duplication_p1 = 5.
 - **normal**: The value will be sampled from a normal distribution between the corresponding values of p0 and p1. For example, if you want the duplication rate is randomly sampled from an normal distribution N(3,5) you should use duplication_d = normal; duplication_p0 = 3; duplication_p1 = 5.

**STEM_FAMILIES**
Number of gene families already in the stem clade

**STEM_LENGTH**
If different from 0, a stem of length STEM_LENGTH is added to the root. In this stem the origination of new gene families occur with a probability proportional to its length. NOT WRITTEN YET!!

**TRANSFER_CLOSE_SPECIES**  

If True (1), transfers take place preferentially between closely related species. When a transfer event occurs, the recipient is chosen randomly among all coexisting lineages with a probability inversely proportional to the logarithm of the evolutionary distance between the two clades.  NOT WRITTEN YET!!

#### Sequence evolution ####

**TL**

**SEQUENCE_EVOLUTION_MODE**

 - 0: No heterogeneity
 - 1: Lineage-specific heterogeneity (time and lineage autocorrelated)
 - 2: Lineage-specific heterogeneity (lineage autocorrelated)
 - 3: Lineage-specific heterogeneity (uncorrelated)
 - 4: User defined heterogeneity

**SUBSTITUTION_D**
**SUBSTITUTION_P0**
**SUBSTITUTION_P1**

