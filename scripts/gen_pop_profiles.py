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
    prob_pass: float = 0.5,
):
    population = [ [int(P(prob_pass)) for _ in range(test_cases)] for _ in range(pop_size) ]
    return population

# TODO - move existing logic here
def GenRandomWithSpecialistsPop(
    pop_size: int,
    num_test_cases: int,
    pass_prob: float,
    ko_prob: float,
    num_specialist_cases: int,
    specialists_per_case: int
):
    num_normal_test_cases = num_test_cases - num_specialist_cases
    num_specialist_candidates = num_specialist_cases * specialists_per_case

    if num_specialist_cases >= num_test_cases:
        print("Cannot have more specialist test cases than total test cases")
        exit(-1)

    if pop_size < num_specialist_candidates:
        print("Pop size too small for configured number of specialists")
        exit(-1)

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

    return population


# Evenly divide population among mutually exclusive niches
def GenEvenExclusiveNichesPop(
    pop_size: int,
    num_test_cases: int,
    num_niches: int,
    pass_prob: float,
    ko_prob: float
):

    if pop_size < num_niches:
        print("Pop size too small for configured number of niches")
        exit(-1)

    if num_test_cases < num_niches:
        print(f"num_test_cases ({num_test_cases}) too small for num_niches ({num_niches})")
        exit(-1)

    # print(f"num_test_cases={num_test_cases}")
    # print(f"num_niches={num_niches}")
    # How many tests per niche?
    tests_per_niche = num_test_cases // num_niches
    num_niche_tests = {}
    unallocated_tests = num_test_cases - (tests_per_niche * num_niches)
    # How many individuals per niche?
    indiv_per_niche = pop_size // num_niches
    num_niche_indivs = {}
    unallocated_indivs = pop_size - (indiv_per_niche * num_niches)
    # Allocate tests and individuals to each niche
    for niche_i in range(num_niches):
        num_niche_tests[niche_i] = tests_per_niche
        num_niche_indivs[niche_i] = indiv_per_niche
        if unallocated_tests > 0:
            num_niche_tests[niche_i] += 1
            unallocated_tests -= 1
        if unallocated_indivs > 0:
            num_niche_indivs[niche_i] += 1
            unallocated_indivs -= 1
    # For each niche, assign individuals with appropriate test case profile
    # print(num_niche_tests)
    start_test = 0
    start_indiv = 0
    population = [None for _ in range(pop_size)]
    for niche_i in range(num_niches):
        # print(f"---niche {niche_i}---")
        # Compute the test profile for this niche
        num_tests = num_niche_tests[niche_i]
        base_test_profile = [i >= start_test and i < (start_test + num_tests) for i in range(num_test_cases)]
        start_test += num_tests
        # Assign individuals this test profile.
        num_indivs = num_niche_indivs[niche_i]
        for indiv_i in range(start_indiv, start_indiv + num_indivs):
            # print(f"  -- indiv {indiv_i} --")
            # Modify base task profile according to pass_prob and ko_prob.
            population[indiv_i] = [int(test and P(pass_prob) and (not P(ko_prob))) for test in base_test_profile]
            # print(f"  profile={population[indiv_i]}")
        start_indiv += num_indivs

    return population

'''
fixed_pass determines if we guarantee a % of passes
pass_prob here determines probability of passing the single discrim case
ko_prob determines whether individuals get background value or 0
'''
def GenSingleDiscrimPop(
    pop_size: int,
    num_test_cases: int,
    fixed_pass: bool,
    pass_prob: float,
    ko_prob: float,
    background: float = 1.0
):
    num_normal_test_cases = num_test_cases - 1
    if fixed_pass:
        num_passes = int(pass_prob * pop_size)
        population = [[background if not P(ko_prob) else 0 for _ in range(num_normal_test_cases)] + [0.0] for _ in range(pop_size)]
        for i in range(num_passes):
            population[i%pop_size][-1] = 1.0
    else:
        population = [[background if not P(ko_prob) else 0 for _ in range(num_normal_test_cases)] + [float(P(pass_prob))] for _ in range(pop_size)]

    return population

# TODO
# def GenCriticalCasesPop():
#     pass


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
    specialists_per_case: int,
    num_niches: int,
    background_value: float,
    fixed_pass: bool
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
            "num_niches": num_niches,
            "background_value": background_value,
            "fixed_pass": fixed_pass,
            "pop_id": pop_id
        }
        population = []
        if pop_type == "random":
            population = GenRandomWithSpecialistsPop(
                pop_size,
                num_test_cases,
                pass_prob,
                ko_prob,
                num_specialist_cases,
                specialists_per_case
            )
        elif pop_type == "even-exclusive-niches":
            population=GenEvenExclusiveNichesPop(
                pop_size,
                num_test_cases,
                num_niches,
                pass_prob,
                ko_prob
            )
        elif pop_type == "single-discrim":
            population = GenSingleDiscrimPop(
                pop_size,
                num_test_cases,
                fixed_pass,
                pass_prob,
                ko_prob,
                background_value
            )
        # print(population)
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
    parser.add_argument("--type", type=str, choices=["random", "even-exclusive-niches", "single-discrim"], help="Type of population profiles to generate")
    parser.add_argument("--pop_size", type=int, default=1000, help="Population size")
    parser.add_argument("--test_cases", type=int, default=100, help="Number of test cases")
    parser.add_argument("--replicates", type=int, default=1, help="")
    parser.add_argument("--pass_prob", type=float, default=0.5, help="Probability candidate passes test case that is randomly assigned as pass/fail."),
    parser.add_argument("--ko_prob", type=float, default=0.0, help="Probability of knocking out a pass once passes have been assigned")
    parser.add_argument("--num_specialist_cases", type=int, default=0, help="Number of test cases that are only solved by specialist specializing on that case.")
    parser.add_argument("--specialists_per_case", type=int, default=1, help="Number of each type of specialist (determined by num_specialist_cases)")
    parser.add_argument("--num_niches", type=int, default=1, help="Number of niches.")
    parser.add_argument("--background_value", type=float, default=1.0, help="Background value to assign tests (single-descrim)" )
    parser.add_argument("--fixed_pass", type=bool, default=True, help="Fix number of discrim test passes")

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
    num_niches = args.num_niches
    fixed_pass = args.fixed_pass
    background_value = args.background_value

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
        specialists_per_case,
        num_niches,
        background_value,
        fixed_pass
    )


if __name__ == "__main__":
    main()