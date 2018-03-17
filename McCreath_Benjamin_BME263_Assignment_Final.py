"""Loads GTF and PSL data and plots them using matplotlib as per BME263 final instructions."""
# links to external data files:
# https://drive.google.com/file/d/1nfuV-67CRPWjUwmr5h-l_D_khIlXlVpS/view?usp=sharing
# https://drive.google.com/file/d/1bZsl9URPTDs_1lIwylw04YlcPRfIfzKl/view?usp=sharing
# https://drive.google.com/file/d/1E-5Opukhtcx0mPQahi49Kggn1788YV_F/view?usp=sharing
import csv
import pickle
import sys
import time
from itertools import groupby
from operator import itemgetter

import matplotlib.patches as mplpatches
import matplotlib.pyplot as plt


def load_data() -> (list, list, list):
    """Load the data from the three specified files."""
    # If the values are passed in via command line, assign them
    if len(sys.argv) == 4:
        gtf_file = sys.argv[1]
        data5 = sys.argv[2]
        data6 = sys.argv[3]
    # If they're not, just hard code them.
    else:
        gtf_file = "gencode.vM12.annotation.gtf"
        data5 = "BME163_Input_data5.psl"
        data6 = "BME163_Input_data6.psl"
    # Initialize the containers
    gtf_data = []
    data5_data = []
    data6_data = []
    # Combine the files with the containers.
    for i, file in enumerate(
            [list(data) for data in zip([gtf_file, data5, data6], [gtf_data, data5_data, data6_data])]):
        # Open files and read them in.
        with open(file[0])as file_in:
            reader = csv.reader(file_in, delimiter="\t")
            if not i:
                file[1] += list(reader)[5:]
            else:
                file[1] += list(reader)
    return gtf_data, data5_data, data6_data


def check_pickle() -> (list, list, list):
    """Check to see if pickled data exists, if it does not, create it. Return all data from pickle or fresh read."""
    # Assuming grader will run off C: drive. I'm running on D: and want pkl. Running on not D: will not give pickle.
    if "D:" in sys.argv[0][:2]:
        try:
            with open("gtf.pkl", mode='r+b') as open_gtf_pickle, \
                    open("data5.pkl", mode='r+b') as open_data5_pickle, \
                    open("data6.pkl", mode='r+b') as open_data6_pickle:
                gtf_data = pickle.load(open_gtf_pickle)
                data5_data = pickle.load(open_data5_pickle)
                data6_data = pickle.load(open_data6_pickle)
        except FileNotFoundError as _:
            gtf_data, data5_data, data6_data = load_data()
            with open("gtf.pkl", mode='w+b') as open_gtf_pickle, \
                    open("data5.pkl", mode='w+b') as open_data5_pickle, \
                    open("data6.pkl", mode='w+b') as open_data6_pickle:
                pickle.dump(gtf_data, open_gtf_pickle)
                pickle.dump(data5_data, open_data5_pickle)
                pickle.dump(data6_data, open_data6_pickle)
    else:
        gtf_data, data5_data, data6_data = load_data()
    return gtf_data, data5_data, data6_data


def process_gtf_data(gtf_data: list) -> list:
    """Extract useful information from gtf data file."""
    clean_data = []
    for i, row in enumerate(gtf_data):
        row[3:5] = list(map(int, row[3:5]))
        # make sure the data is in the desired range
        if 'chr7' in row and any([45232945 <= row[3] <= 45240000, 45232945 <= row[4] <= 45240000]):
            if any([row[2] == "transcript", row[2] == "exon", row[2] == "CDS"]):
                row[-1:] = row[-1].split(";")
                name = [x.strip() for x in row if "transcript_name" in str(x)][0].split('\"')[1]
                clean_row = [name] + list(itemgetter(2, 3, 4)(row))
                # clean_row[0] = float(clean_row[0].rsplit('G')[1].replace("\"", ""))
                clean_data.append(clean_row)
    clean_data.sort(key=itemgetter(3), reverse=True)
    # sortdict captures the order of the transcript IDs when sorted by end index so the transcript IDs can be sorted
    # by index instead of string value
    sortdict = {}
    i = 0
    for row in clean_data:
        if row[0] not in sortdict:
            sortdict[row[0]] = i
            i += 1
    # See? It's sorted down here.
    clean_data = sorted(clean_data, key=lambda x: sortdict[x[0]], reverse=True)
    return clean_data


def process_psl(raw_psl_data: list) -> list:
    """Extract useful data from raw PSL data."""
    psl_data = []
    for row in raw_psl_data:
        clean_row = list(itemgetter(13, 15, 16, 20, 18)(row))
        clean_row[1:3] = list(map(int, clean_row[1:3]))
        # Make sure it's in the right area and just strip the data.
        if 'chr7' in clean_row and any([45232945 <= clean_row[1] <= 45240000, 45232945 <= clean_row[2] <= 45240000]):
            for num in [3, 4]:
                clean_row[num] = list(map(int, clean_row[num].split(',')[:-1]))
            psl_data.append(clean_row)
    psl_data.sort(key=itemgetter(2))
    return psl_data


def process_data(gtf_data: list, data5_data: list, data6_data: list) -> (list, list, list):
    """Drive the processing of the data in one function call."""
    # Apply each data file to the correct processing function.
    return process_gtf_data(gtf_data), process_psl(data5_data), process_psl(data6_data)


def plot_gtf(panel: plt.Axes, data: list, y_offset: float, linewidths: list) -> plt.Axes:
    """Apply gtf data to specified data axis."""
    # initialize y_index
    y_index = y_offset
    # This dictionary is used to plot hte width based on type
    lwdd = {"transcript": linewidths[0], "exon": linewidths[1], "CDS": linewidths[2]}
    # This dictionary holds the last valued plotted at each y_index.
    y_index_dict = {}
    # Groups by each transcript ID
    for _, group in groupby(data, itemgetter(0)):
        # Convert from iterator to list.
        group = list(group)
        new_start = min([x[2] for x in group])
        new_end = max([x[3] for x in group])
        # set the y_index
        if y_index_dict.keys():
            y_index = max(y_index_dict.keys()) + y_offset
        for dict_key in y_index_dict:
            if new_start > y_index_dict[dict_key]:
                y_index = min(dict_key, y_index)
        for feature in group:
            panel.add_patch(mplpatches.Rectangle(xy=(feature[2], (y_index - (lwdd[feature[1]] / 2))),
                                                 width=(feature[3] - feature[2]),
                                                 height=lwdd[feature[1]], edgecolor=None, facecolor='Black'))
        y_index_dict[y_index] = new_end
    return panel


def plot_psl(panel: plt.Axes, data: list, y_offset: float, linewidths: list) -> plt.Axes:
    """Apply psl data to specified data axis."""
    # initialize the offset
    y_index = y_offset
    # Dict holds Y values and furthest right point of any line.
    y_index_dict = {}
    # checks to see if line overlaps other line(s) and sets y index accordingly
    for _, record in enumerate(data):
        if y_index_dict.keys():
            y_index = max(y_index_dict.keys()) + y_offset
        for dict_key in y_index_dict:
            if record[1] > y_index_dict[dict_key]:
                y_index = min(dict_key, y_index)
        # Adds the main rectangle
        panel.add_patch(mplpatches.Rectangle(xy=(record[1], y_index), width=record[2] - record[1],
                                             height=linewidths[0], edgecolor=None, facecolor='Black'))
        # Adds all the features of the main rectangle.
        for feature in zip(record[3], record[4]):
            panel.add_patch(mplpatches.Rectangle(xy=(feature[0], (y_index - (linewidths[1] / 2))), width=feature[1],
                                                 height=linewidths[1], edgecolor=None, facecolor='Black'))
        # updates dictionary
        y_index_dict[y_index] = record[2]
    return panel


def main():
    """Drive plot creation and set properties. Load data and such as necessary."""
    start = time.time()
    gtf_data, data5_data, data6_data = process_data(*check_pickle())
    plt.style.use('BME163.mplstyle')
    plt_size = [10, 5]
    plt.figure(figsize=plt_size)
    # Create the three panels.
    top = plt.axes([0 / plt_size[0], 3.25 / plt_size[1], 10 / plt_size[0], 1.25 / plt_size[1]], xticks=[], yticks=[],
                   xlim=(45232945, 45240001), ylim=(0, 1))
    mid = plt.axes([0 / plt_size[0], 1.75 / plt_size[1], 10 / plt_size[0], 1.25 / plt_size[1]], xticks=[], yticks=[],
                   xlim=(45232945, 45240001), ylim=(0, 1))
    bot = plt.axes([0 / plt_size[0], 0.25 / plt_size[1], 10 / plt_size[0], 1.25 / plt_size[1]], xticks=[], yticks=[],
                   xlim=(45232945, 45240001), ylim=(0, 1))
    # Plot the three data files.
    plot_gtf(top, gtf_data, .1, [.011, .025, .0501])
    plot_psl(mid, data6_data, 1 / 67, [.001, .0075])
    plot_psl(bot, data5_data, .0025, [.0001, .0007])
    plt.savefig('McCreath_Benjamin_BME263_Assignment_Final.png', dpi=1200)
    print("Time to complete: {}s".format(time.time() - start))
    return 0


if __name__ == '__main__':
    main()
