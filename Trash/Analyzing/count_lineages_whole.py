import sys

def count_lineages_in_time_whole(events_file):

    lineages_counter = 0
    output = list()

    with open(events_file) as f:

        f.readline()

        for line in f:

            time, e, nodes = line.strip().split("\t")

            if e == "S":
                lineages_counter += 1

            if e == "E":
                lineages_counter -=1

            if e == "F":

                final_time = float(time)
                break

            output.append((time, str(lineages_counter)))

    for item in output:
        line = "\t".join(item)
        print(line)



if __name__ == "__main__":

    my_file = sys.argv[1:][0]
    count_lineages_in_time_whole(my_file)
