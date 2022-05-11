import argparse

def CalcLexicaseEvals_NoSampling(generations, pop_size, num_tests):
    return generations * pop_size * num_tests

def CalcLexicaseEvals_RandomSampling(generations, pop_size, num_tests, test_ds_rate):
    evals_per_gen = pop_size * int((num_tests * test_ds_rate))
    return generations * evals_per_gen

def CalcLexicaseEvals_MaxMinPopSampling(
    generations,
    pop_size,
    num_tests,
    test_ds_rate,
    pop_sample_rate
):
    # calculate size of population sample
    pop_sample_size = int(pop_size*pop_sample_rate)
    # calculate size of test sample
    test_sample_size = int(num_tests * test_ds_rate)
    # calculate evaluations per generation as number of evals for fully evaluating population sample plus
    # evaluating the rest of the population on the sampled test cases
    evals_per_gen = (pop_sample_size*num_tests) + ((pop_size-pop_sample_size)*test_sample_size)
    return generations*evals_per_gen


def FullLexGensToRandDSLexGens(
    full_lex_gens,
    ds_rate
):
    return full_lex_gens / ds_rate

"""
Convert number of generations of full lexicase to number of generations for
maxmin-pop-sampling down-sampling lexicase.
"""
def FullLexGensToMaxMinPopDSLexGens(
    full_lex_gens,
    pop_size,
    num_tests,
    test_ds_rate,
    pop_sample_rate
):
    numerator = CalcLexicaseEvals_NoSampling(full_lex_gens, pop_size, num_tests)
    pop_sample_size = int(pop_size * pop_sample_rate)
    test_sample_size = int(num_tests * test_ds_rate)
    denom = (pop_sample_size*num_tests) + (test_sample_size * (pop_size-pop_sample_size))
    return numerator / denom



def main():
    parser = argparse.ArgumentParser(description="Calculate number of generations to run for different lexicase sampling regimes.")
    parser.add_argument("--pop_size", type=int, help="Population size")
    parser.add_argument("--num_tests", type=int, help="Number of test cases")
    parser.add_argument("--full_lex_gens", type=int, help="Number of full lexicase generations")
    parser.add_argument("--test_sample_rate", type=float, default=0.5, help="Test case sampling rate (between 0 and 1)")
    parser.add_argument("--pop_sample_rate", type=float, default=0.5, help="Population sampling rate (for maxmin test sampling")

    args = parser.parse_args()
    pop_size = args.pop_size
    num_tests = args.num_tests
    full_lex_gens = args.full_lex_gens
    test_sample_rate = args.test_sample_rate
    pop_sample_rate = args.pop_sample_rate

    full_lex_evals = CalcLexicaseEvals_NoSampling(full_lex_gens, pop_size, num_tests)

    random_sampling_gens = FullLexGensToRandDSLexGens(full_lex_gens, 0.2)
    random_sampling_evals = CalcLexicaseEvals_RandomSampling(random_sampling_gens, pop_size, num_tests, test_sample_rate)

    maxmin_sampling_gens = FullLexGensToMaxMinPopDSLexGens(
        full_lex_gens=full_lex_gens,
        pop_size=pop_size,
        num_tests=num_tests,
        test_ds_rate=test_sample_rate,
        pop_sample_rate=pop_sample_rate
    )
    maxmin_sampling_evals = CalcLexicaseEvals_MaxMinPopSampling(
        maxmin_sampling_gens,
        pop_size,
        num_tests,
        test_sample_rate,
        pop_sample_rate
    )

    print("Full lexicase:")
    print(f"  Generations: {full_lex_gens}")
    print(f"  Evaluations: {full_lex_evals}")

    print("Random down-sampling lexicase:")
    print(f"  Generations: {random_sampling_gens}")
    print(f"  Evaluations: {random_sampling_evals}")

    print("Maxmin down-sampling lexicase:")
    print(f"  Generations: {maxmin_sampling_gens}")
    print(f"  Evaluations: {maxmin_sampling_evals}")

if __name__ == "__main__":
    main()