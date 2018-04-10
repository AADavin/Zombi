from SpeciesTreeSimulator import SpeciesTreeGenerator
from GenomeSimulator import GenomeSimulator
from SequenceSimulator import SequenceSimulator
import AuxiliarFunctions as af
import argparse
import os
import sys

class simuLyon():

    def __init__(self):

        self.tree_parameters = dict()
        self.genome_parameters = dict()
        self.sequence_parameters = dict()


    def T(self, parameters_file, experiment_folder, advanced_mode):

        tree_folder = os.path.join(experiment_folder, "T")

        if advanced_mode == "i":

            # In this case the input is a tree file
            tree_file = parameters_file
            print("Generate events for input file %s" % tree_file)

            stg = SpeciesTreeGenerator({})

            stg.start()

            stg.events = af.generate_events(tree_file)

            events_file = os.path.join(tree_folder, "Events.tsv")
            stg.write_events_file(events_file)

            whole_tree_file = os.path.join(tree_folder, "WholeTree.nwk")
            extant_tree_file = os.path.join(tree_folder, "ExtantTree.nwk")

            whole_tree, extant_tree = stg.generate_newick_trees()

            with open(whole_tree_file, "w") as f:
                f.write(whole_tree)

            with open(extant_tree_file, "w") as f:
                f.write(extant_tree)

            return 0

        else:

            parameters = af.prepare_species_tree_parameters(af.read_parameters(parameters_file))
            stg = SpeciesTreeGenerator(parameters)

        os.system("cp " + parameters_file + " " + tree_folder)

        run_counter = 0
        success = False

        while success == False and run_counter <= 50:
            run_counter+=1
            print("Computing Species Tree. Trial number %s" % str(run_counter))
            if advanced_mode == "0":
                success = stg.run()
            if advanced_mode == "a":
                success = stg.run_a()
            if advanced_mode == "b":
                success = stg.run_b()
            if advanced_mode == "p":
                success = stg.run_p()

        if run_counter >= 100:
            print("Aborting computation of the Species Tree. Please use other speciation and extinction rates!")
            return 0

        events_file = os.path.join(tree_folder, "Events.tsv")
        stg.write_events_file(events_file)

        whole_tree_file = os.path.join(tree_folder, "WholeTree.nwk")
        extant_tree_file = os.path.join(tree_folder, "ExtantTree.nwk")

        whole_tree, extant_tree = stg.generate_newick_trees()

        with open(whole_tree_file, "w") as f:
            f.write(whole_tree)

        with open(extant_tree_file, "w") as f:
            f.write(extant_tree)

        if advanced_mode == "b" or advanced_mode == "a":
            rates_file = os.path.join(tree_folder, "Rates.tsv")
            stg.write_rates(rates_file)


    def G(self,parameters_file, experiment_folder, advanced_mode):

        parameters = af.prepare_genome_parameters(af.read_parameters(parameters_file))
        events_file = os.path.join(experiment_folder, "T/Events.tsv")
        genome_folder = os.path.join(experiment_folder, "G")
        os.system("cp " + parameters_file + " " + genome_folder)

        genomes_folder = os.path.join(genome_folder, "Genomes")
        gene_families_folder = os.path.join(genome_folder, "Gene_families")

        gss = GenomeSimulator(parameters, events_file)

        if advanced_mode == "0":
            gss.run()

        elif advanced_mode == "u":

            rates_folder = os.path.join(experiment_folder, "CustomRates")
            gss.read_rates(rates_folder)
            gss.run_u()

        # We write the output

        print("Writing Genomes")
        gss.write_genomes(genomes_folder)
        print("Writing Gene Families")
        gss.write_gene_family_events(gene_families_folder)

        if parameters["PROFILES"] == 1:
            print("Writing Profiles")
            profiles_folder = os.path.join(genome_folder, "Profiles")
            gss.write_profiles(profiles_folder)

        if parameters["EVENTS_PER_BRANCH"] == 1:
            print("Writing Events Per Branch")
            events_per_branch_folder = os.path.join(genome_folder, "Events_per_branch")
            gss.write_events_per_branch(events_per_branch_folder)

        if parameters["GENE_TREES"] == 1:
            print("Writing Gene Trees")
            gene_trees_folder = os.path.join(genome_folder, "Gene_trees")
            gss.write_gene_trees(gene_trees_folder)

    def S(self, parameters_file, experiment_folder, advanced_mode):

        gene_trees_folder = os.path.join(experiment_folder, "G/Gene_trees")
        sequences_folder = os.path.join(experiment_folder, "S/Sequences")
        os.system("cp " + parameters_file + " " + sequences_folder)
        fasta_folder = os.path.join(sequence_folder, "Sequences")

        if not os.path.isdir(sequences_folder):
            os.mkdir(sequences_folder)

        if not os.path.isdir(fasta_folder):
            os.mkdir(fasta_folder)

        parameters = af.prepare_sequence_parameters(af.read_parameters(parameters_file))

        print("Preparing simulator of sequences")

        ss = SequenceSimulator(parameters)

        if advanced_mode == "0":

            whole_trees = [x.replace("_pruned", "_whole") for x in os.listdir(gene_trees_folder) if "pruned" in x]

            for tree_file in whole_trees:
                tree_path = os.path.join(gene_trees_folder, tree_file)
                print("Simulating sequence for gene family %s" % tree_file.split("_")[0])
                ss.run(tree_path, fasta_folder)
                ss.write_pruned_sequences(tree_path.replace("whole", "pruned"), fasta_folder)

        elif advanced_mode == "u":

            # First we obtain the rates-multiplier

            ss.obtain_rates_multipliers(experiment_folder + "/CustomRates/GT_Substitution_rates.tsv",
                                        experiment_folder + "/CustomRates/ST_Substitution_rates.tsv")

            # And we save it

            ss.write_rates_sttree(experiment_folder + "/T/WholeTree.nwk",
                                  os.path.join(experiment_folder, "T/RatesTree.nwk"))

            whole_trees = [x for x in os.listdir(gene_trees_folder) if "whole" in x]


            for tree_file in whole_trees:
                tree_path = os.path.join(gene_trees_folder, tree_file)
                print("Simulating sequence for gene family %s" % tree_file.split("_")[0])
                ss.run_u(tree_path, fasta_folder)
                ss.write_pruned_sequences(tree_path.replace("whole", "pruned"), fasta_folder)












if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("mode", type=str, choices=["T","Ti","Tb","Tp","G", "Gu","S","Su"], help="Mode")
    parser.add_argument("params",  type=str, help="Parameters file")
    parser.add_argument("output", type=str, help="Name of the experiment folder")

    args = parser.parse_args()

    mode, parameters_file, experiment_folder =  args.mode, args.params, args.output

    if len(mode) == 1:
        main_mode = mode[0]
        advanced_mode = "0"

    elif len(mode) == 2:
        main_mode = mode[0]
        advanced_mode = mode[1]

    else:
        print ("Incorrect value for mode")

    SL = simuLyon()

    if main_mode == "T":

        if not os.path.isdir(experiment_folder):
            os.mkdir(experiment_folder)

        if not os.path.isdir(os.path.join(experiment_folder,"T")):
            os.mkdir(os.path.join(experiment_folder, "T"))

        SL.T(parameters_file, experiment_folder, advanced_mode)

    elif main_mode == "G":

        genome_folder = os.path.join(experiment_folder, "G")

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)

        SL.G(parameters_file, experiment_folder, advanced_mode)

    elif main_mode == "S":

        sequence_folder = os.path.join(experiment_folder, "S")

        if not os.path.isdir(sequence_folder):
            os.mkdir(sequence_folder)

        SL.S(parameters_file, experiment_folder, advanced_mode)


    else:
        print("Incorrect usage. Please select a mode: T, G or S")

