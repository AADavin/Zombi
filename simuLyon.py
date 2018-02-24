from GeneFamilyGenerator import FamilyOriginator, GeneFamilySimulator
from GenomeGenerator import GenomeSimulator
from TreeGenerator import TreeGenerator
from TreeGeneratorContinuousTime import TreeGeneratorContinuousTime
import ete3
import sys
import os

class SimuLYON():

    def __init__(self):

        self.tree_parameters = dict()
        self.genome_parameters = dict()
        self.sequence_parameters = dict()

    def read_parameters(self, parameters_file, parameters):

        with open(parameters_file) as f:
            for line in f:
                if line[0] == "#":
                    continue
                parameter, value = line.strip().split("\t")
                parameters[parameter] = value

    def obtain_species_trees(self, parameters_file, experiment_folder):

        self.read_parameters(parameters_file, self.tree_parameters)

        if self.tree_parameters["SPECIES_EVOLUTION_MODE"] == '0':

            success = False
            trials = 0

            while success == False and trials <= 100:
                print("Computing Species Tree, %s trial" % str(trials))
                tg = TreeGenerator(parameters_file)
                success = tg.new_tree_generator()
                trials += 1

            if success == False:
                print("Aborting. Maximum number of trials attained. Not possible to compute the tree")

            else:
                print("Correctly computed tree")
                tg.store_log(experiment_folder)

        # You have to work in the next lines!

        elif self.tree_parameters["SPECIES_EVOLUTION_MODE"] == '1':
            tg.generate_tree_mode_1()
        elif self.tree_parameters["SPECIES_EVOLUTION_MODE"] == '2':
            tg.generate_tree_mode_2()
        elif self.tree_parameters["SPECIES_EVOLUTION_MODE"] == '3':
            tg.generate_tree_mode_3()
        elif self.tree_parameters["SPECIES_EVOLUTION_MODE"] == '4':
            tg.generate_tree_mode_4()

    def obtain_species_trees_continuous_time(self, parameters_file, experiment_folder):

        self.read_parameters(parameters_file, self.tree_parameters)

        if self.tree_parameters["SPECIES_EVOLUTION_MODE"] == '0':

            success = False
            trials = 0

            while success == False and trials <= 100:
                print("Computing Species Tree, %s trial" % str(trials))
                tg = TreeGeneratorContinuousTime(parameters_file)
                success = tg.new_tree_generator()
                trials += 1

            if success == False:
                print("Aborting. Maximum number of trials attained. Not possible to compute the tree")

            else:
                print("Correctly computed tree")
                tg.store_log(experiment_folder)


    def obtain_gene_families(self, parameters_file, experiment_folder):

        self.read_parameters(parameters_file, self.genome_parameters)

        prefix = self.genome_parameters["PREFIX"]

        stopping_rule = int(self.genome_parameters["STOPPING_RULE"])
        families_in_stem = int(self.genome_parameters["STEM_FAMILIES"])
        n_families = int(self.genome_parameters["N_FAMILIES"])

        whole_tree_file = os.path.join(experiment_folder, "WholeTree")

        with open(whole_tree_file) as f:
            tree = ete3.Tree(f.readline().strip(),format=1)

        events_file = os.path.join(experiment_folder, "SpeciesTreeEvents.tsv")
        lineages_file = os.path.join(experiment_folder, "LineagesInTime.tsv")

        gene_families_folder = os.path.join(experiment_folder, "GeneFamilies")

        if not os.path.isdir(gene_families_folder):
           os.mkdir(gene_families_folder)

        profiles_file = os.path.join(gene_families_folder, "Profiles.tsv")
        transfers_file = os.path.join(gene_families_folder, "Transfers.tsv")

        with open(transfers_file, "w") as f:
            f.write("Family\tDonor\tRecipient\n")

        self._prepare_profile_file(tree, profiles_file)

        fo = FamilyOriginator(whole_tree_file, events_file)
        gfs = GeneFamilySimulator(parameters_file, events_file, lineages_file)

        if stopping_rule == 0:

            for i in range(1, families_in_stem + 1):

                family_name = prefix + str(i)

                print("Simulating family %s" % family_name)
                gfs.origination("Root", 0, family_name)
                gfs.run_mode_0()
                gfs.complete_gene_family_information()

                tree = gfs.get_gene_family_tree()
                if tree != "None":
                    with open(os.path.join(gene_families_folder, family_name), "w") as f:
                        f.write(tree)

                gfs.output_profile(profiles_file, "Profile")
                gfs.write_transfers(transfers_file)
                gfs.write_log(os.path.join(gene_families_folder, family_name + "_events.tsv"))

            for j in range(i+1, n_families + i + 1):

                family_name = prefix + str(j)

                node, time = fo.create_families()

                print("Simulating family %s appearing in node %s and in time %s" % (family_name, node, str(time)))

                gfs.origination(node, time, family_name)
                gfs.run_mode_0()
                gfs.complete_gene_family_information()

                tree = gfs.get_gene_family_tree()

                if tree != "None":
                    with open(os.path.join(gene_families_folder, family_name), "w") as f:
                        f.write(tree)

                gfs.output_profile(profiles_file, "Profile")
                gfs.write_transfers(transfers_file)
                gfs.write_log(os.path.join(gene_families_folder, family_name + "_events.tsv"))


    def obtain_genomes(self, parameters_file, experiment_folder):

        self.read_parameters(parameters_file, self.genome_parameters)

        whole_tree_file = os.path.join(experiment_folder, "WholeTree")
        events_file = os.path.join(experiment_folder, "SpeciesTreeEvents.tsv")
        lineages_file = os.path.join(experiment_folder, "LineagesInTime.tsv")
        profiles_file = os.path.join(experiment_folder, "Profiles.tsv")
        transfers_file = os.path.join(experiment_folder, "Transfers.tsv")

        with open(transfers_file, "w") as f:
            f.write("Family\tDonor\tRecipient\n")

        genome_folder = os.path.join(experiment_folder, "Genomes")
        events_per_branch_folder = os.path.join(genome_folder, "EventsPerBranch")
        gene_trees_folder = os.path.join(genome_folder, "GeneTrees")
        events_per_family_folder = os.path.join(genome_folder, "EventsPerFamily")
        complete_genomes_folder = os.path.join(genome_folder, "CompleteGenomes")

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)
        if not os.path.isdir(events_per_branch_folder):
            os.mkdir(events_per_branch_folder)
        if not os.path.isdir(events_per_family_folder):
            os.mkdir(events_per_family_folder)
        if not os.path.isdir(gene_trees_folder):
            os.mkdir(gene_trees_folder)
        if not os.path.isdir(complete_genomes_folder):
            os.mkdir(complete_genomes_folder)

        with open(whole_tree_file) as f:
            tree = ete3.Tree(f.readline().strip(),format=1)

        self._prepare_profile_file(tree, profiles_file)

        fo = FamilyOriginator(whole_tree_file, events_file)
        gfs = GenomeSimulator(parameters_file, events_file, lineages_file)

        gfs.run(complete_genomes_folder )

        print("Writing logs")

        gfs.write_log(gene_trees_folder, events_per_family_folder, complete_genomes_folder, events_per_branch_folder)


    def obtain_sequences(self, parameters_file, experiment_folder):
        pass

    def _prepare_profile_file(self,tree, myfile):

        with open(myfile, "w") as f:
            line = list()
            line.append("Family")
            myroot = tree.get_tree_root()
            myroot.name = "Root"
            for node in tree.traverse():
                line.append(str(node.name))
            myline = "\t".join(line) + "\n"
            f.write(myline)

if __name__ == "__main__":

    args = sys.argv[1:]
    if len(args) != 3:
        print("Incorrect usage. Please read the manual. The usual way to run this script is:")
        print("python simuLyon.py T Parameters_file.tsv /Output_folder")
    else:
        mode, parameters_file, experiment_folder = args

        SL = SimuLYON()

        if mode == "T":

            if os.path.isdir(experiment_folder):
                pass
            else:
                os.mkdir(experiment_folder)
            SL.obtain_species_trees(parameters_file, experiment_folder)

        if mode == "TC":

            if os.path.isdir(experiment_folder):
                pass
            else:
                os.mkdir(experiment_folder)
            SL.obtain_species_trees_continuous_time(parameters_file, experiment_folder)

        elif mode == "G":

            SL.obtain_genomes(parameters_file, experiment_folder)

        elif mode == "F":

            SL.obtain_gene_families(parameters_file, experiment_folder)

        elif mode == "S":

            print("This mode is not ready to be used yet")

        else:
            print("Incorrect usage. Please select a mode: T, G, F or S")
