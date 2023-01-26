# 2023-01-26-even-niche-static-analysis experiment

Goal: analyze ability of different sampling methods to select initial wave of individuals that solve a new test case.

General parameters:

- Pop Size: 1000
- Tests: 200
- Selection: random down-sampling, informed down-sampling, standard lexicase
- DS rate: 0.05, 0.1, 0.2
- Pop sampling (informed): 0.01, 0.1

Synthetic population parameters

- Type: single-discrim
- BG = 0.0
- pass_prob = 0.01, 0.1, 0.2
- fixed_pass= true


20 replicates of each.
