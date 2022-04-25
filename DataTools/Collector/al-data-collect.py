#####################################################################################################
#####################################################################################################
# Create csv's with requested data from all replicates every X generations
#
#
# Output : csv with for data over time
#
# python3
#####################################################################################################
#####################################################################################################

######################## IMPORTS ########################
import pandas as pd
import numpy as np
import argparse
import math as mth
import sys
import os

# variables we are testing for each replicate range
TR_LIST = ['1','2','4','8','16','32','64','128','256']
TS_LIST = ['1','2','4','8','16','32','64','128','256']
LX_LIST = ['0.0','0.1','0.3','0.6','1.2','2.5','5.0']
FS_LIST = ['0.0','0.1','0.3','0.6','1.2','2.5','5.0']
ND_LIST = ['0.0','0.1','0.3','0.6','1.2','2.5','5.0']
NS_LIST = ['0','1','2','4','8','15','30']

# seed experiements replicates range
REP_NUM = 50

# columns we are interested in grabbing
# pop level
POP_FIT_AVG = 'pop_fit_avg'
POP_FIT_MAX = 'pop_fit_max'
POP_OPT_AVG = 'pop_opt_avg'
POP_OPT_MAX = 'pop_opt_max'
POP_UNI_OBJ = 'pop_uni_obj'
POP_STR_AVG = 'pop_str_avg'
POP_STR_MAX = 'pop_str_max'
# elite solutions (max agg traits)
ELE_AGG_PER = 'ele_agg_per'
ELE_OPT_CNT = 'ele_opt_cnt'
# # optimal solution (max opti trait count)
# OPT_AGG_PER = 'opt_agg_per'
# OPT_OBJ_CNT = 'opt_obj_cnt'
# # streak solution (longest streak of non-zero values)
# STR_AGG_PER = 'str_agg_per'
# STR_OBJ_CNT = 'str_obj_cnt'
# trait coverage
UNI_STR_POS = 'uni_str_pos'
# # pareto data
PARETO_CNT = 'pareto_cnt'
# # novelty data
ARCHIVE_CNT = 'archive_cnt'
PMIN = 'pmin'
# ARC_ELITE = 'arc_elite'
# ARC_OPTI = 'arc_opti'
# ARC_STRK = 'arc_strk'
ARC_ACTI_GENE = 'arc_acti_gene'
OVERLAP = 'overlap'

# return appropiate string dir name (based off run.sb file naming system)
def SetSelection(s,p):
    # case by case
    if s == 0:
        return 'TRUNCATION'
    elif s == 1:
        return 'TOURNAMENT'
    elif s == 2:
        if p == '0':
            return 'FITSHARING_G'
        elif p == '1':
            return 'FITSHARING_P'
        else:
            sys.exit("UNKNOWN SELECTION")
    elif s == 3:
        return 'LEXICASE'
    elif s == 4:
        return 'NONDOMINATEDSORTING'
    elif s == 5:
        return 'NOVELTY'
    else:
        sys.exit("UNKNOWN SELECTION")

# return appropiate string dir name (based off run.sb file naming system)
def SetDiagnostic(s):
    # case by case
    if s == 0:
        return 'EXPLOITATION'
    elif s == 1:
        return 'STRUCTEXPLOITATION'
    elif s == 2:
        return 'STRONGECOLOGY'
    elif s == 3:
        return 'EXPLORATION'
    elif s == 4:
        return 'WEAKECOLOGY'
    else:
        sys.exit('UNKNOWN DIAGNOSTIC')

# return appropiate string dir name (based off run.sb file naming system)
def SetSelectionVar(s):
    # case by case
    if s == 0:
        return 'TR'
    elif s == 1:
        return 'T'
    elif s == 2:
        return 'SIG'
    elif s == 3:
        return 'EPS'
    elif s == 4:
        return 'SIG'
    elif s == 5:
        return 'NOV'
    else:
        sys.exit("UNKNOWN SELECTION VAR")

# return the correct amount of seed ran by experiment treatment
def SetSeeds(s):
    # case by case
    if s == 0 or s == 1:
        seed = []
        seed.append([x for x in range(1,51)])
        seed.append([x for x in range(51,101)])
        seed.append([x for x in range(101,151)])
        seed.append([x for x in range(151,201)])
        seed.append([x for x in range(201,251)])
        seed.append([x for x in range(251,301)])
        seed.append([x for x in range(301,351)])
        seed.append([x for x in range(351,401)])
        seed.append([x for x in range(401,451)])
        return seed

    elif s == 2 or s == 3 or s == 4 or s == 5:
        seed = []
        seed.append([x for x in range(1,51)])
        seed.append([x for x in range(51,101)])
        seed.append([x for x in range(101,151)])
        seed.append([x for x in range(151,201)])
        seed.append([x for x in range(201,251)])
        seed.append([x for x in range(251,301)])
        seed.append([x for x in range(301,351)])
        return seed

    else:
        sys.exit('SEEDS SELECTION UNKNOWN')

# set the appropiate list of variables we are checking for
def SetVarList(s):
    # case by case
    if s == 0:
        return TR_LIST
    elif s == 1:
        return TS_LIST
    elif s == 2:
        return FS_LIST
    elif s == 3:
        return LX_LIST
    elif s == 4:
        return FS_LIST
    elif s == 5:
        return NS_LIST
    else:
        sys.exit("UNKNOWN VARIABLE LIST")

# return extra parameter directory if needed
def SetSecondParam(s, pt):
    # case by case
    if s == 0:
        return ''
    elif s == 1:
        return ''
    elif s == 2:
        return ''
    elif s == 3:
        return ''
    elif s == 4:
        return ''
    elif s == 5:
        return ''
    else:
        sys.exit("UNKNOWN SELECTION")

# loop through differnt files that exist
def DirExplore(data, dump, sel, dia, offs, res, obj, acc, gens, pt):
    # check if data dir exists
    if os.path.isdir(data) == False:
        print('DATA=', data)
        sys.exit('DATA DIRECTORY DOES NOT EXIST')

    # check if data dir exists
    if os.path.isdir(dump) == False:
        print('DATA=', data)
        sys.exit('DATA DIRECTORY DOES NOT EXIST')

    # check that selection data folder exists
    SEL_DIR = data + SetSelection(sel,pt) + '/TRT_' + obj + '__ACC_' + acc + '__GEN_' + gens + '/'
    if os.path.isdir(SEL_DIR) == False:
        print('SEL_DIR=', SEL_DIR)
        sys.exit('SELECTION DIRECTORY NOT FOUND')

    # Set vars that we need to loop through
    VLIST = SetVarList(sel)
    SEEDS = SetSeeds(sel)

    # gens we are expecting
    GEN_LIST = [x for x in range(int(gens)+1) if x%res == 0]
    # data frame list for concatanation
    DF_LIST = []
    # second parameter dir
    SECOND_PARAM = SetSecondParam(sel, pt)

    # iterate through the sets of seeds
    for i in range(len(SEEDS)):
        seeds = SEEDS[i]
        print('i=',i)

        TRT = [VLIST[i]] * len(GEN_LIST)

        # iterate through seeds to collect data
        for s in seeds:
            it = int((s-1)/REP_NUM)
            var_val = str(VLIST[it])
            seed = str(s + offs)
            DATA_DIR =  SEL_DIR + 'DIA_' + SetDiagnostic(dia) + '__' + SetSelectionVar(sel) + '_' + var_val + '__SEED_' + seed + '/' + SECOND_PARAM + 'data.csv'

            # check if data file even exists
            if os.path.isfile(DATA_DIR) == False:
                print('DNE: ' + DATA_DIR)
                continue

            # create pandas data frame of entire csv and grab the row
            df = pd.read_csv(DATA_DIR)
            df = df.iloc[::res, :]

            # time to export the data
            cdf = pd.DataFrame({'gen': pd.Series(GEN_LIST),
                            'trt': pd.Series(TRT),
                            'pop_fit_avg': pd.Series(df[POP_FIT_AVG].tolist()),
                            'pop_fit_max': pd.Series(df[POP_FIT_MAX].tolist()),
                            'pop_opt_avg': pd.Series(df[POP_OPT_AVG].tolist()),
                            'pop_opt_max': pd.Series(df[POP_OPT_MAX].tolist()),
                            'pop_uni_obj': pd.Series(df[POP_UNI_OBJ].tolist()),
                            'pop_str_avg': pd.Series(df[POP_STR_AVG].tolist()),
                            'pop_str_max': pd.Series(df[POP_STR_MAX].tolist()),
                            'ele_agg_per': pd.Series(df[ELE_AGG_PER].tolist()),
                            'ele_opt_cnt': pd.Series(df[ELE_OPT_CNT].tolist()),
                            'arc_acti_gene': pd.Series(df[ARC_ACTI_GENE].tolist()),
                            # 'opt_agg_per': pd.Series(df[OPT_AGG_PER].tolist()),
                            # 'opt_obj_cnt': pd.Series(df[OPT_OBJ_CNT].tolist()),
                            # 'str_agg_per': pd.Series(df[STR_AGG_PER].tolist()),
                            # 'str_obj_cnt': pd.Series(df[STR_OBJ_CNT].tolist()),
                            'overlap': pd.Series(df[OVERLAP].tolist()),
                            'archive_cnt': pd.Series(df[ARCHIVE_CNT].tolist()),
                            'uni_str_pos':  pd.Series(df[UNI_STR_POS].tolist()),
                            'pmin':        pd.Series(df[PMIN].tolist()),
                            'pareto_cnt' :  pd.Series(df[PARETO_CNT].tolist())}) #),
                            # 'arc_elite':   pd.Series(df[ARC_ELITE].tolist()),
                            # 'arc_opti':    pd.Series(df[ARC_OPTI].tolist()),
                            # 'arc_strk':    pd.Series(df[ARC_STRK].tolist())})

            DF_LIST.append(cdf)

    fin_df = pd.concat(DF_LIST)

    fin_df.to_csv(path_or_buf= dump + 'al-' + SetDiagnostic(dia).lower() + '-' + gens + '-' + obj + '-' + acc + '.csv', index=False)

def main():
    # Generate and get the arguments
    parser = argparse.ArgumentParser(description="Data aggregation script.")
    parser.add_argument("data_dir",    type=str, help="Target experiment directory.")
    parser.add_argument("dump_dir",    type=str, help="Data dumping directory")
    parser.add_argument("selection",      type=int, help="Selection scheme we are looking for? \n0: Truncation\n1: Tournament\n2: Fitness Sharing\n4: Espilon Lexicase\n6: Nondominated Sorting\n7: Novelty Search")
    parser.add_argument("diagnostic",     type=int, help="Diagnostic we are looking for?\n0: Exploitation\n1: Structured Exploitation\n2: Strong Ecology\n3: Exploration\n4: Weak Ecology")
    parser.add_argument("seed_offset", type=int, help="Experiment seed offset. (REPLICATION_OFFSET + PROBLEM_SEED_OFFSET")
    parser.add_argument("resolution",  type=int, help="The resolution desired for the data extraction")
    parser.add_argument("objectives", type=str, help="Number of objectives being optimized")
    parser.add_argument("accuracy", type=str, help="Accuracy for experiment")
    parser.add_argument("generations", type=str, help="Number of generations experiments ran for")
    parser.add_argument("param_two",      type=str, help="Second paramater for any selection scheme")

    # Parse all the arguments
    args = parser.parse_args()
    data_dir = args.data_dir.strip()
    print('Data directory=',data_dir)
    dump_dir = args.dump_dir.strip()
    print('Dump directory=', dump_dir)
    selection = args.selection
    print('Selection scheme=', SetSelection(selection, args.param_two))
    diagnostic = args.diagnostic
    print('Diagnostic=',SetDiagnostic(diagnostic))
    offset = args.seed_offset
    print('Offset=', offset)
    resolution = args.resolution
    print('Resolution=', resolution)
    objectives = args.objectives
    print('Objectives=', objectives)
    accuracy = args.accuracy
    print('Accuracy=', accuracy)
    generations = args.generations
    print('Generations=', generations)
    param_two = args.param_two
    print('2nd param=', param_two)

    # Get to work!
    print("\nChecking all related data directories now!")
    DirExplore(data_dir, dump_dir, selection, diagnostic, offset, resolution, objectives, accuracy, generations, param_two)

if __name__ == "__main__":
    main()