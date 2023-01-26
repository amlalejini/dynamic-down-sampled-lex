#!/usr/bin/env bash

REPLICATES=20
EXP_SLUG=2023-01-26-specialist-static-analysis
ACCOUNT=devolab
SEED_OFFSET=5000
JOB_TIME=24:00:00
JOB_MEM=16G

SCRATCH_EXP_DIR=./test/data/dynamic-ds-lex
REPO_DIR=/home/amlalejini/devo_ws/dynamic-down-sampled-lex
HOME_EXP_DIR=${REPO_DIR}/experiments

DATA_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}
JOB_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}/jobs
CONFIG_DIR=${HOME_EXP_DIR}/${EXP_SLUG}/hpc/config

python3 gen-sub.py --time_request ${JOB_TIME} --mem ${JOB_MEM} --data_dir ${DATA_DIR} --config_dir ${CONFIG_DIR} --repo_dir ${REPO_DIR} --replicates ${REPLICATES} --job_dir ${JOB_DIR} --account ${ACCOUNT} --seed_offset ${SEED_OFFSET}