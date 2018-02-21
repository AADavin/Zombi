
##  
##### A new tool to simulate genome evolution accounting for dead lineages

----------

### **Introduction** ###

There are two methods to simulate the evolution of genomes

The first one is simulating the evolution of gene families independently, one by one. Genes families evolve inside the species tree undergoing events 
of duplications, transfers and losses. Simulating the evolution of many gene families one can simulate the evolution of genomes. 
This method however, does not provide any information regarding the position of the different genes within the genomes. Genomes are simply a bunch
of genes that are found in no particular order in a given point of the species tree. This method however, has the advantage that you can use **gene-wise** rates (easier to estimate) and that you can simulate as many gene families as you want

The second method starts with an ancestral genome at the root, with a given number of genes. For now in the current version, all genes families present in this 
genome have a single copy (so in this ancestral genome there are no duplicated genes). In the current version, genes have 4 fields, separated by underscored.
For example:

n4_+_5_12

The meaning of this is:
n4: The branch in the species tree where that gene is found
+: The orientation of the gene (it can be + or -). This is determined randomly with a bernouilli distribution p=0.5
5: The gene family identifier (in this case, gene family 5)
12: An unique identifier for that gene. 
 
A genome is an ordered collection of genes. So if we begin with a genome that has 5 genes, what we see is something like

Position    Gene_family Orientation Id
0           1           +           1
1           2           -           1
2           3           -           1
3           4           +           1
4           5           -           1

The genome is circular, so the position 4 is adjacent to the position 3 and 0

The events of evolution in this mode are:

D: Duplications. A segment of the genome is duplicated. The new copy is inserted next to the old one
L: Losses. A segment of the genome is lost
T: Transfers. A segment of the genome is transferred to a contemporary species. The segment is inserted in a random position
C: Translocations. A segment of the genome changes its position within the genome
I: Inversions. A segment of the genome inverts its position
O: Originations. A new gene family appears and it is inserted in a random position

The rates in this case are **genome-wise**.
There is also an additional rate for each event. This is called the **extension_rate**. This number (between 0 and 1) is
the p parameters of a **geometric distribution** that controls the length of the affected segment.

For example, if DUPLICATION_E == 1, the extension of the segment duplicated will be always 1 (meaning that only one gene is duplicated at a time)

By changing this parameter we can fine tune the extension associated to the different events. If inversions affect normally large chunks of the genome,
it suffices to use a low p.

Origination of new gene families are always of size 1, meaning that it is not possible to have an origination of two gene families in the same step of time.

Once that the full evolution of genomes has been simulated, simuLyon prints also the gene trees associated to the different gene families, all the events taking place
in each gene family, the events taking place in each branch and the genomes of each node in the species tree.










