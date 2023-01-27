'''
Aggregate data
'''

import argparse, os, sys, pathlib
sys.path.append(os.path.join(pathlib.Path(os.path.dirname(os.path.abspath(__file__))).parents[2], "scripts"))
import utilities as utils

run_identifier = "RUN_"


cfg_fields_exclude = [
    "OUTPUT_DIR"
]
data_fields_exclude = [
    "pop_info",
    "pop_test_coverage_profile"
]

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")

    # Parse command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    dump_dir = args.dump

    # Verify that the given data directory exits
    if not os.path.exists(data_dir):
        print("Unable to find data directory.")
        exit(-1)

    # Create the directory to dump aggregated data (if it doesn't already exist)
    utils.mkdir_p(dump_dir)

    # Aggregate run directories.
    run_dirs = [run_dir for run_dir in os.listdir(data_dir) if run_identifier in run_dir]
    print(f"Found {len(run_dirs)} run directories.")

    # Create file to hold time series data
    time_series_header = None   # Holds the time series file header (verified for consistency across runs)
    time_series_fpath = os.path.join(dump_dir, f"time_series.csv")
    with open(time_series_fpath, "w") as fp:
        fp.write("")

    incomplete_runs = []
    total_runs = len(run_dirs)
    cur_run_i = 0
    for run_dir in run_dirs:
        run_path = os.path.join(data_dir, run_dir)

        cur_run_i += 1
        print(f"Processing ({cur_run_i}/{total_runs}) {run_path}")

        # Skip over (but make note of) incomplete runs.
        required_files = [
            os.path.join("output", "summary.csv"),
            os.path.join("output", "run_config.csv")
        ]
        incomplete = any([not os.path.exists(os.path.join(run_path, req)) for req in required_files])
        if incomplete:
            print("  - Failed to find all required files!")
            incomplete_runs.append(run_dir)
            continue

        ############################################################
        # Extract run parameters (from run_config.csv)
        cfg_path = os.path.join(run_path, "output", "run_config.csv")
        cfg_data = utils.read_csv(cfg_path)
        run_cfg = {line["parameter"]:line["value"] for line in cfg_data}
        ############################################################

        ############################################################
        # Extract run data
        run_data_path =  os.path.join(run_path, "output", "summary.csv")
        run_data = utils.read_csv(run_data_path)
        run_info = []
        for line in run_data:
            # Copy cfg details into info
            info = {field:run_cfg[field] for field in run_cfg if field not in cfg_fields_exclude}
            # Copy line data into info
            for field in line:
                if field in data_fields_exclude: continue
                info[field] = line[field]
            # Compose custom fields
            info["pop_class"] = info["POPULATIONS_FILE"].replace(".csv", "")
            info["pop_test_coverage_profile"] = f'\"{line["pop_test_coverage_profile"]}\"'
            # - Build selection condition field
            sel_condition = info["SELECTION_METHOD"]
            test_sample_method = info["TEST_SAMPLING_METHOD"]
            test_sample_prop = info["TEST_SAMPLING_PROP"]
            if test_sample_method!="none":
                sel_condition += f"_{test_sample_method}_t{test_sample_prop}"
            if test_sample_method == "maxmin-sample":
                sel_condition += f"_pop{info['MAXMIN_POP_PROP']}"
            if info["PARTITIONING_METHOD"] != "none":
                sel_condition += f"_{info['PARTITIONING_METHOD']}"
            if info["PARTITIONING_METHOD"] =="random-cohort":
                sel_condition += f"_{info['COHORT_PARTITIONING_PROP']}"
            info["sel_condition"] = sel_condition
            run_info.append(info)

            sampling_method = test_sample_method
            if "random-cohort" in sel_condition:
                sampling_method = "cohort-partitioning"
            sample_prop = info["COHORT_PARTITIONING_PROP"] if sampling_method == "cohort-partitioning" else info["TEST_SAMPLING_PROP"]

            info["sample_method"] = sampling_method
            info["sample_prop"] = sample_prop
        ############################################################

        ############################################################
        # Write output
        # Analyze fields, build header if necessary
        fields = list(run_info[0].keys())
        fields.sort()
        write_header = False
        if time_series_header==None:
            time_series_header = fields
            write_header = True
        elif time_series_header != fields:
            print("Header mismatch!")
            exit(-1)
        out_content = "\n".join([",".join([info[key] for key in fields]) for info in run_info])
        with open(time_series_fpath, "a") as fp:
            if write_header:
                fp.write(",".join(time_series_header) + "\n")
            fp.write(out_content)
            fp.write("\n")
        ############################################################



if __name__ == "__main__":
    main()
