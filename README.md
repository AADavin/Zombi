

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

[Link to the paper](https://watermark.silverchair.com/btz710.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAq8wggKrBgkqhkiG9w0BBwagggKcMIICmAIBADCCApEGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMeuxJkjvPpGqRHF73AgEQgIICYmMCnTbF25wuitWpQMCfPQmPqS9VBWumzI0oQ8Cd8qf-lZ9AO1BfOp6AphbqD9HGqJHniR6gDnjxze2O2t0pLReW9mJ27Zev2OVx5wLluGysTxK3ddTN7UT81omAUQj43IGkOMF3mE039d-dwleAcT2OX0iXAR8kZoOGRXIl13PYEiPeZqai1SUbGfd9plkECzJbhFaHKb-SvOAPzyxwaNDhDrRPWxaTdAAnhcSLL3w-uUHwW43P91gHXUHrpDxd1K5U-sGZx7Ej538eJ6tVKhO_GEyLuyFqpNFKnOyAvtw_364ivK6DWyFblvLI4yS275Iq3g77EZ8Ezm7NJxdn4-2tiy84X_uUKijMewA3AxWdOyA7u8WqBYQ6khBDeiN9rw5v8QaCNGc4ZWy41psqxqKobXx4fIjrvNBXBT05XF6IQ_YkyMXcwfgVMewmeO_VsQftgJ4_yO3fTvr3qZlqbcTzLs_6W9Kj3WU8xcK_MuXQr_hUaI6DnFY5MG6aAXtS8cc6_DyzAddHIkZ7qroF52m-hIEz4GhBU_Cou16ftDygxJQbUQMLHFeLPQv1dJFWZVyDwQvj14-3tmZpf8jq7tOs4PrqnrUrcB6X6ZpnHEdSEL-T1eJeRoHuDIlbm8pHhvxRB0DSuOjktR54h0R1fX7uppgRdLvC6tu8Qu6-KG6q23nzrjXQ2XvI8z5EikM46QL4RQaZN58yVuKowXMB-4EdaNsePr89XQ37Q9yc18o9uiSHqVAiry5-T3FEd_2EhmqnsEraHkLvduchvBp4-w94BlY9utjTGFYs-n0cyOw4xGg)


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


