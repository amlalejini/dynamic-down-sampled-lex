#!/bin/bash

SCRIPT=../../../scripts/gen_pop_profiles.py
SEED=2
POP_SIZE=1000
TEST_CASES=200
NUM_SPECIALIST_CASES=0
SPECIALISTS_PER_CASE=0
PASS_PROB=1.0
KO_PROB=0.0
BG_VALUE=0.0
FIXED_PASS=1

TYPE=even-exclusive-niches

NUM_NICHES=10
OUTPUT=even-excl-niches-${NUM_NICHES}.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}

NUM_NICHES=20
OUTPUT=even-excl-niches-${NUM_NICHES}.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}

NUM_NICHES=50
OUTPUT=even-excl-niches-${NUM_NICHES}.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}