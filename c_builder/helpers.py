"""Helper functions for the parsing of files/directories."""
# Directory finding/parsing
from os import listdir
from os.path import join
from os.path import isdir

# Regex parsing
from re import match


def pick_sample(file_struct):
    """Choose a sample directory from the file_struct."""
    for base_key in file_struct:
        if base_key == "name":
            continue
        for middle_key in file_struct[base_key]:
            if middle_key == "name":
                continue
            for end_key in file_struct[base_key][middle_key]:
                if end_key == "name":
                    continue
                return file_struct[base_key][middle_key][end_key]


def parse_directory(dir_path):
    """Recursively parse the contents of a directory."""
    dir_struct = {}

    # Build Directory structure
    for file_ in listdir(dir_path):
        file_path = join(dir_path, file_)
        if isdir(file_path):
            # Find structure of subdirectory
            dir_struct[file_] = parse_directory(file_path)
            dir_struct[file_]["name"] = file_
        elif file_path.endswith(('.mhd', '.nii')):
            # We only care about .mhd and .nii files
            dir_struct[file_] = {}
            dir_struct[file_]['file'] = 1
            dir_struct[file_]["name"] = file_

    return dir_struct


def parse_c(c_csv, skip_header=True):
    """Parse contents of C file and return dict."""
    # Skip header row
    if skip_header:
        next(c_csv)

    data = {}
    # Go through csv file
    for row in c_csv:
        if row[0] == ' ' or row[0] == '':
            # Empty rows
            continue

        # Get values for file
        incremental = row[3]
        modality = row[4]
        label = row[6]
        reg = row[7]

        # File path parsing
        file_path = row[8]
        file_path_parts = file_path.split("/")
        first_folder = file_path_parts[0]
        middle_folder = file_path_parts[1]
        last_folder = file_path_parts[2]
        filename = "/".join(file_path_parts[3:])

        # Build resulting dictionary with the file info
        if first_folder not in data:
            data[first_folder] = {}
        if middle_folder not in data[first_folder]:
            data[first_folder][middle_folder] = {}
        if last_folder not in data[first_folder][middle_folder]:
            data[first_folder][middle_folder][last_folder] = {}

        # Modality is ins/exp
        # Reg is F01/L01F01
        # Label is filename/data/R

        # Set the values for this file
        this_datum = {}
        this_datum["file"] = 1
        this_datum["incremental"] = incremental
        this_datum["mod"] = modality
        this_datum["label"] = label
        this_datum["reg"] = reg
        this_datum["name"] = filename

        data[first_folder][middle_folder][last_folder][filename] = this_datum

    return data


def generate_csv_rows(data_dict, sample_data, csv_out, first, middle, last, subdir_name=""):
    """Given an object write the data to csv."""
    name = "/".join([subdir_name, data_dict["name"]])
    if data_dict.get("file", 0) == 1:  # We have a file
        # Default row
        out_row = [first, middle, last, "FIX ME!",
                   "FIX ME!", "FIX ME!", "FIX ME!", "FIX ME!", "FIX ME!"]

        for sample_file in sample_data:
            sample_file_data = sample_data[sample_file]
            if match(sample_file_data["regex"], name):
                if sample_file_data["ignore"]:
                    return

                if not sample_file_data["incremental"]:
                    out_row[3] = 0
                else:
                    out_row[3] = 1

                if not sample_file_data["modality"] == "":
                    out_row[4] = sample_file_data["modality"]

                out_row[5] = subdir_name

                if not sample_file_data["label"] == "":
                    out_row[6] = sample_file_data["label"]

                if not sample_file_data["reg"] == "":
                    out_row[7] = sample_file_data["reg"]

                break

        out_row[8] = "{}/{}/{}{}".format(first, middle, last, name)
        csv_out.writerow(out_row)

    else:  # We have a directory
        for file_info in data_dict:
            if file_info == "name":
                continue

            generate_csv_rows(data_dict[file_info],
                              sample_data, csv_out, first, middle, last, subdir_name=name)
