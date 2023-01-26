#!/bin/bash

SCRIPT=../../../scripts/gen_pop_profiles.py
SEED=2
POP_SIZE=1000
TEST_CASES=200
KO_PROB=0.0
BG_VALUE=1.0
FIXED_PASS=1
NUM_NICHES=1
PASS_PROB=1.0

TYPE=random
NUM_SPECIALIST_CASES=1

SPECIALISTS_PER_CASE=200
OUTPUT=specialists-t${NUM_SPECIALIST_CASES}-p${SPECIALISTS_PER_CASE}.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}

SPECIALISTS_PER_CASE=100
OUTPUT=specialists-t${NUM_SPECIALIST_CASES}-p${SPECIALISTS_PER_CASE}.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}

SPECIALISTS_PER_CASE=50
OUTPUT=specialists-t${NUM_SPECIALIST_CASES}-p${SPECIALISTS_PER_CASE}.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}

SPECIALISTS_PER_CASE=10
OUTPUT=specialists-t${NUM_SPECIALIST_CASES}-p${SPECIALISTS_PER_CASE}.csv
python3 ${SCRIPT} --dump ./dump/ --output ${OUTPUT} --seed ${SEED} --type ${TYPE} --pop_size ${POP_SIZE} --test_cases ${TEST_CASES} --replicates 1 --pass_prob ${PASS_PROB} --ko_prob ${KO_PROB} --num_specialist_cases ${NUM_SPECIALIST_CASES} --specialists_per_case ${SPECIALISTS_PER_CASE} --num_niches ${NUM_NICHES} --background_value ${BG_VALUE} --fixed_pass ${FIXED_PASS}