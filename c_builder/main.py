"""Server file for the program."""
# CSV parsing
from csv import reader
from csv import writer

# Request data handling
import json

# File loading methods
from Tkinter import Tk
from tkFileDialog import askdirectory
from tkFileDialog import askopenfile
from tkFileDialog import asksaveasfilename

# Flask imports to make it work as a website
from flask import Flask
from flask import render_template
from flask import jsonify
from flask import request

# Helper functions
from helpers import pick_sample
from helpers import parse_c
from helpers import parse_directory
from helpers import generate_csv_rows
from helpers import sort_keys

C_BUILDER = Flask(__name__)

DATA = None
timepoints = None


@C_BUILDER.route('/')
def index():
    """Return the base page."""
    global DATA
    DATA = None
    return render_template('index.html')


@C_BUILDER.route('/files/')
def load_files():
    """Handle request to get a directory."""
    print("Loading files")
    global DATA, timepoints

    timepoints = request.args.get("timepoint") == "true"

    # Hide default window that opens
    root = Tk()
    # root.withdraw()
    root.update()

    # Prompt user for directory
    folder = askdirectory()

    # Clean up tkinter window
    root.destroy()
    
    # Cancel button clicked, abort
    if not folder:
        print("Load Aborted")
        return jsonify({'status': False})

    # Directory selected, get folder structure
    output = {}
    folder_struct = parse_directory(folder)
    if timepoints:
        folder_struct = {folder: folder_struct}

    DATA = folder_struct

    output['sample'] = pick_sample(folder_struct)
    output['status'] = True

    # Return results
    return jsonify(output)


@C_BUILDER.route('/c/')
def load_c():
    """Handle request to parse C file."""
    global DATA

    # Hide default window
    root = Tk()
    # root.withdraw()
    root.update()

    c_file = askopenfile()
    root.destroy()
    if c_file is None:
        print("Load Aborted")
        return jsonify({'status': False})

    c_csv = reader(c_file)

    output = {}
    parsed_c = parse_c(c_csv)
    output['status'] = True
    output['sample'] = pick_sample(parsed_c)
    DATA = parsed_c

    return jsonify(output)


@C_BUILDER.route('/export/', methods=['POST'])
def export_file():
    """Export data to csv file."""
    print("Exporting data")
    global DATA

    if DATA is None:
        print("Directory not set.")
        return jsonify({"status": False})

    data_in = request.form.to_dict(flat=False)
    data_in = list(data_in.keys())[0]

    sample_data = json.loads(data_in)

    # Hide default window that opens
    root = Tk()
    # root.withdraw()
    root.update()

    out_filename = asksaveasfilename()
    root.destroy()

    if not out_filename:
        print("Export aborted!")
        return jsonify({'status': False})

    if "." not in out_filename:
        out_filename = "{}.csv".format(out_filename)

    # Weird fix for windows outputting multiple newlines
    csv_out = writer(open(out_filename, mode='w'))
    csv_out.writerow(["first folder", "middle", "last folder", "Incremental[1]",
                      "Modality[2]", "subdirectory", "data/label/R/jac", "Reg", "Location"])

    for first_folder in sort_keys(DATA):
        if first_folder == "name":
            continue
        first_folder_data = DATA[first_folder]
        for second_folder in sort_keys(first_folder_data):
            if second_folder == "name":
                continue
            second_folder_data = first_folder_data[second_folder]
            for last_folder in sort_keys(second_folder_data):
                if last_folder == "name":
                    continue
                last_folder_data = second_folder_data[last_folder]
                for file_ in sort_keys(last_folder_data):
                    if file_ == "name":
                        continue
                    file_data = last_folder_data[file_]
                    first_folder_name = first_folder.split("/")[-1]
                    generate_csv_rows(file_data, sample_data, csv_out, first_folder_name, second_folder, last_folder)

    response = {}
    response["status"] = True
    return jsonify(response)


@C_BUILDER.route('/clear/')
def clear_file():
    """Clear stored data."""
    print("Clearing data")
    global DATA
    DATA = None
    return jsonify({"status": True})
