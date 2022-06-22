# 2022-06-21-sel-analysis

Run preliminary version of SelectorAnalyzer to compare different test case subsampling approaches on synthetic populations with randomly generated test case profiles (pass/fail assignments).

- Selection methods:
  - Lexicase
  - Lexicase + random test sample
  - Lexicaxe + maxmin test sample (1%, 10%)
  - Lexicase + cohort partitioning
  - Tournament
  - Tournament + random test sample
  - Tournament + maxmin test sample (1%, 10%)
  - Tournament + cohort partitioning
- Sample rates:
  - 5%, 10%
- Gens
  - 300
- Networks
  - N=10, Pop=1000, tests=200, type=Random, pass probability=50%, specialists=10
  - N=10, Pop=1000, tests=200, type=Random, pass probability=10%, specialists=10
  - N=10, Pop=1000, tests=200, type=Random, pass probability=5%, specialists=10
  - N=10, Pop=1000, tests=200, type=Random, pass probability=1%, specialists=10

Commands used to generate networks:

```
python gen_pop_profiles.py --dump ./pop_profiles/ --type random --pop_size 1000 --test_cases 200 --replicates 10 --pass_prob 0.5 --ko_prob 0.0 --num_specialist_cases 10 --specialists_per_case 1 --output random-1000p-200t-50pp-10s.csv
python gen_pop_profiles.py --dump ./pop_profiles/ --type random --pop_size 1000 --test_cases 200 --replicates 10 --pass_prob 0.1 --ko_prob 0.0 --num_specialist_cases 10 --specialists_per_case 1 --output random-1000p-200t-10pp-10s.csv
python gen_pop_profiles.py --dump ./pop_profiles/ --type random --pop_size 1000 --test_cases 200 --replicates 10 --pass_prob 0.05 --ko_prob 0.0 --num_specialist_cases 10 --specialists_per_case 1 --output random-1000p-200t-05pp-10s.csv
python gen_pop_profiles.py --dump ./pop_profiles/ --type random --pop_size 1000 --test_cases 200 --replicates 10 --pass_prob 0.01 --ko_prob 0.0 --num_specialist_cases 10 --specialists_per_case 1 --output random-1000p-200t-01pp-10s.csv
```
