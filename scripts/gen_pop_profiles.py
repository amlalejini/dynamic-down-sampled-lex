"""
Generate populations of pass/fail test case performance profiles.

The generated test case performance profiles can be used to analyze differences in which candidates selection schemes choose.
"""

import argparse, os, random, pyvarco
import utilities

def P(p):
    return random.random() < p

def GenRandomPop(
    pop_size: int,
    test_cases: int,
    prob_pass: float = 0.5
):
    population = [ [int(P(prob_pass)) for _ in range(test_cases)] for _ in range(pop_size) ]
    return population

def gen_pop_profiles(
    dump_dir: str,
    output: str,
    pop_type: str,
    pop_size: int,
    num_test_cases: int,
    replicates: int,
    pass_prob: float,
    ko_prob: float,
    num_specialist_cases: int,
    specialists_per_case: int
):
    output_fname = output
    # Setup output file
    if output_fname == None:
        output_fname = f"type-{pop_type}.csv"
    utilities.mkdir_p(dump_dir)
    output_fpath = os.path.join(dump_dir, output_fname)

    output_content = []
    output_header = None
    write_header = True
    pop_id_set = set()

    # If output file exists, parse header.
    if os.path.exists(output_fpath):
        with open(output_fpath, "r") as fp:
            lines = fp.read().split("\n")
            output_header = lines[0].split(",")
        # Populate pop_id_set
        prev_content = utilities.read_csv(output_fpath)
        pop_id_set = {int(line["pop_id"]) for line in prev_content}
        write_header = False
    else:
        # Otherwise, create the output file
        with open(output_fpath, "w") as fp:
            fp.write("")

    num_normal_test_cases = num_test_cases - num_specialist_cases
    num_specialist_candidates = num_specialist_cases * specialists_per_case

    if num_specialist_cases >= num_test_cases:
        print("Cannot have more specialist test cases than total test cases")
        exit(-1)

    if pop_size < num_specialist_candidates:
        print("Pop size too small for configured number of specialists")
        exit(-1)

    prev_pop_id = -1 if len(pop_id_set) == 0 else max(pop_id_set)
    for rep_i in range(replicates):
        pop_id = (prev_pop_id + 1)
        prev_pop_id += 1

        rep_info = {
            "pop_type": pop_type,
            "pop_size": pop_size,
            "num_test_cases": num_test_cases,
            "pass_prob": pass_prob,
            "ko_prob": ko_prob,
            "num_specialist_cases": num_specialist_cases,
            "specialists_per_case": specialists_per_case,
            "pop_id": pop_id
        }

        if pop_type == "random":
            population = GenRandomPop(
                pop_size = pop_size-num_specialist_candidates,
                test_cases = num_normal_test_cases,
                prob_pass=pass_prob
            )

        # Run probabilistic knockouts on each individual's test cases
        if ko_prob > 0.0:
            for cand_i in range(len(population)):
                for test_i in range(len(population[cand_i])):
                    if (population[cand_i][test_i] > 0) and P(ko_prob):
                        population[cand_i][test_i]=0

        # Add specialist test cases to current population
        for cand in population:
            cand += [0 for _ in range(num_specialist_cases)]

        # Add specialist individuals to the population
        # specialist_pop_ids = {i for i in range(len(population), pop_size)}
        for spec_test_i in range(0, num_specialist_cases):
            population += [ [int(test_i == num_normal_test_cases+spec_test_i) for test_i in range(0, num_test_cases)] for cand_i in range(0, specialists_per_case) ]

        ################################################
        # Add population to output content
        for cand_i in range(len(population)):
            info = {field:rep_info[field] for field in rep_info}
            candidate = population[cand_i]
            info["candidate_id"] = cand_i
            info["test_scores"] =  '"[' + ",".join(list(map(str, candidate))) + ']"' #f'"{str(candidate)}"'
            output_content.append(info)
        ################################################

    fields = list(output_content[0].keys())
    fields.sort()
    output_lines = [ ",".join([str(line[field]) for field in fields]) for line in output_content]
    if not write_header and fields != output_header:
        print("Header mismatch:", output_header, fields)

    with open(output_fpath, "a") as fp:
        if write_header:
            fp.write(",".join(fields))
        fp.write("\n")
        fp.write("\n".join(output_lines))

def main():
    parser = argparse.ArgumentParser(description="Generate population test case profiles.")
    parser.add_argument("--dump", type=str, help="Where to dump profiles?", default=".")
    parser.add_argument("--output", type=str, default=None, help="File name of output. If not given, construct a default file name.")
    parser.add_argument("--seed", type=int, default=2, help="Random number seed")
    parser.add_argument("--type", type=str, choices=["random", "niches"], help="Type of population profiles to generate")
    parser.add_argument("--pop_size", type=int, default=1000, help="Population size")
    parser.add_argument("--test_cases", type=int, default=100, help="Number of test cases")
    parser.add_argument("--replicates", type=int, default=1, help="")
    parser.add_argument("--pass_prob", type=float, default=0.5, help="Probability candidate passes test case that is randomly assigned as pass/fail."),
    parser.add_argument("--ko_prob", type=float, default=0.0, help="Probability of knocking out a pass once passes have been assigned")
    parser.add_argument("--num_specialist_cases", type=int, default=0, help="Number of test cases that are only solved by specialist specializing on that case.")
    parser.add_argument("--specialists_per_case", type=int, default=1, help="Number of each type of specialist (determined by num_specialist_cases)")

    # Parse command line arguments
    args = parser.parse_args()
    random_seed = args.seed
    dump_dir = args.dump
    pop_type = args.type
    pop_size = args.pop_size
    num_test_cases = args.test_cases
    pass_prob = args.pass_prob
    ko_prob = args.ko_prob
    num_specialist_cases = args.num_specialist_cases
    specialists_per_case = args.specialists_per_case
    replicates = args.replicates
    output_fname = args.output

    # Seed the random number generator
    random.seed(random_seed)

    gen_pop_profiles(
        dump_dir,
        output_fname,
        pop_type,
        pop_size,
        num_test_cases,
        replicates,
        pass_prob,
        ko_prob,
        num_specialist_cases,
        specialists_per_case
    )

if __name__ == "__main__":
    main()