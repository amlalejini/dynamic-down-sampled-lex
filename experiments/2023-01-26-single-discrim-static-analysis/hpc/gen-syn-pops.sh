#!/bin/bash

SCRIPT=../../../scripts/gen_pop_profiles.py
SEED=2
POP_SIZE=1000
TEST_CASES=200
NUM_SPECIALIST_CASES=0
SPECIALISTS_PER_CASE=0
KO_PROB=0.0
BG_VALUE=0.0
FIXED_PASS=1
NUM_NICHES=1

TYPE=single-discrim

PASS_PROB=0.01
OUTPUT=single-discrim-01.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}

PASS_PROB=0.05
OUTPUT=single-discrim-05.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}

PASS_PROB=0.1
OUTPUT=single-discrim-10.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}

PASS_PROB=0.2
OUTPUT=single-discrim-20.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}