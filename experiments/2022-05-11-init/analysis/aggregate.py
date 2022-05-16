'''
Aggregate data
'''

import argparse, os, sys, pathlib
sys.path.append(os.path.join(pathlib.Path(os.path.dirname(os.path.abspath(__file__))).parents[2], "scripts"))
import utilities as utils

run_identifier = "RUN_"

time_series_data_fields = [
    "gen",
    "evals",
    "pop_uni_obj",
    "ele_agg_per",
    "ele_opt_cnt"
]
time_series_cfg_fields = [
    "POP_SIZE",
    "MAX_GENS",
    "SEED",
    "SELECTION",
    "DIAGNOSTIC",
    "LEX_DS_RATE",
    "LEX_DS_MODE",
    "condition_name"
]

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--units", type=str, default="gen", choices=["gen", "evals", "interval", "total"], help="Unit for resolution of time series")
    parser.add_argument("--resolution", type=int, default=1, help="What resolution should we collect time series data at?")

    # Parse command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    dump_dir = args.dump
    # update = args.update
    update = -1
    # time_series_range = args.time_series_range

    time_series_units = args.units
    time_series_resolution = args.resolution

    # Verify that the given data directory exits
    if not os.path.exists(data_dir):
        print("Unable to find data directory.")
        exit(-1)

    # Verify time series resolution >= 1
    if time_series_resolution < 1:
        print("Time series resolution must be >= 1")
        exit(-1)

    # Create the directory to dump aggregated data (if it doesn't already exist)
    utils.mkdir_p(dump_dir)

    # Aggregate run directories.
    run_dirs = [run_dir for run_dir in os.listdir(data_dir) if run_identifier in run_dir]
    print(f"Found {len(run_dirs)} run directories.")

    # Create file to hold time series data
    time_series_content = []    # This will hold all the lines to write out for a single run; written out for each run.
    time_series_header = None   # Holds the time series file header (verified for consistency across runs)
    time_series_fpath = os.path.join(dump_dir, f"time_series.csv")
    time_series_update_set = None
    with open(time_series_fpath, "w") as fp:
        fp.write("")

    # Create data structure to hold summary information for each run (1 element per run)
    summary_header = None
    summary_content_lines = []

    incomplete_runs = []

    # Loop over runs, aggregating data from each.
    total_runs = len(run_dirs)
    cur_run_i = 0
    for run_dir in run_dirs:
        run_path = os.path.join(data_dir, run_dir)

        summary_info = {}                   # Hold summary information about run. Indexed by field.
        time_series_info = {}               # Hold time series information. Indexed by update.

        cur_run_i += 1
        print(f"Processing ({cur_run_i}/{total_runs}): {run_path}")

        # Skip over (but make note of) incomplete runs.
        required_files = [
            os.path.join("output", "data.csv"),
            os.path.join("output", "elite.csv"),
            os.path.join("cmd.log")
        ]
        incomplete = any([not os.path.exists(os.path.join(run_path, req)) for req in required_files])
        if incomplete:
            print("  - Failed to find all required files!")
            incomplete_runs.append(run_dir)
            continue

        ############################################################
        # Extract commandline configuration settings (from cmd.log file)
        cmd_log_path = os.path.join(run_path, "cmd.log")
        cmd_params = utils.extract_params_cmd_log(cmd_log_path)
        print(cmd_params)

        condition_name = f'{cmd_params["SELECTION"]}_s-{cmd_params["LEX_DS_MODE"]}_r-{cmd_params["LEX_DS_RATE"]}'
        summary_info["condition_name"] = condition_name
        for field in cmd_params:
            summary_info[field] = cmd_params[field]
        ############################################################

        ############################################################
        # Extract run data
        run_data_path = os.path.join(run_path, "output", "data.csv")
        run_data = utils.read_csv(run_data_path)

        # Identify final generation
        generations = [int(line["gen"]) for line in run_data]
        final_gen = max(generations)

        # Isolate data from final generation
        final_data = [line for line in run_data if int(line["gen"]) == final_gen]
        assert len(final_data) == 1
        final_data = final_data[0]

        # Add final data to run summary info
        for field in final_data:
            summary_info[field] = final_data[field]


        filtered_ts_data = utils.filter_ordered_data(
            run_data,
            time_series_units,
            time_series_resolution
        )
        ts_generations_included = [int(line["gen"]) for line in filtered_ts_data]
        ts_generations_included.sort()
        for line in filtered_ts_data:
            gen = int(line["gen"])
            time_series_info[gen] = {}
            for field in time_series_data_fields:
                time_series_info[gen][field] = line[field]
            for field in time_series_cfg_fields:
                time_series_info[gen][field] = summary_info[field]
        ############################################################

        ############################################################
        # Output time series data for this run
        # Compute time series header from time_series_info
        time_series_fields = list(time_series_info[ts_generations_included[0]].keys())
        time_series_fields.sort()
        # If we haven't written the header, write it.
        write_header = False
        if time_series_header == None:
            write_header = True
            time_series_header = ",".join(time_series_fields)
        elif time_series_header != ",".join(time_series_fields):
            print("Time series header mismatch!")
            exit(-1)

        # Write time series content line-by-line
        time_series_content = []
        for u in ts_generations_included:
            time_series_content.append(",".join([str(time_series_info[u][field]) for field in time_series_fields]))
        with open(time_series_fpath, "a") as fp:
            if write_header: fp.write(time_series_header)
            fp.write("\n")
            fp.write("\n".join(time_series_content))
        time_series_content = []
        ############################################################

        ############################################################
        # Add summary_info to aggregate content
        summary_fields = list(summary_info.keys())
        summary_fields.sort()
        if summary_header == None:
            summary_header = summary_fields
        elif summary_header != summary_fields:
            print("Header mismatch!")
            exit(-1)
        summary_line = [str(summary_info[field]) for field in summary_fields]
        summary_content_lines.append(",".join(summary_line))
        ############################################################

    # write out aggregate data
    with open(os.path.join(dump_dir, "aggregate.csv"), "w") as fp:
        out_content = ",".join(summary_header) + "\n" + "\n".join(summary_content_lines)
        fp.write(out_content)

    # Write out incomplete runs, sort them!
    incomplete_runs.sort()
    with open(os.path.join(dump_dir, "incomplete_runs_agg.log"), "w") as fp:
        fp.write("\n".join(incomplete_runs))

if __name__ == "__main__":
    main()
