import pandas as pd
import os
import csv
import re
import sys
from collections import defaultdict
from math import sqrt, log2, pow

def get_quant_col_prefix(quantification):
    if quantification == "LFQ":
        return "LFQ intensity "
    elif quantification == "intensity":
        return "Intensity "
    elif quantification == "spc":
        return "MS/MS count "

class ProteinGroups:
    def __init__(self, experimental_design, file_path, quantification_saint, quantification_comppass):
        print("\nParsing MaxQuant proteinGroups: " + file_path)

        self.experimental_design = experimental_design
        self.quantification_saint = quantification_saint
        self.quantification_comppass = quantification_comppass
        self.data = pd.read_csv(file_path, sep="\t", quotechar="'", low_memory=False)

        self.quant_col_prefix_saint = get_quant_col_prefix(quantification_saint)
        self.quant_col_prefix_comppass = get_quant_col_prefix(quantification_comppass)

        self.quant_cols_saint = []
        self.quant_cols_comppass = []
        self.bait_cols = []

        # make sure experimental design matches the data
        for col in self.data.columns:
            if col.startswith(self.quant_col_prefix_saint):
                exp_name = col.split(self.quant_col_prefix_saint)[-1]
                if exp_name not in self.experimental_design.name2experiment:
                    print("WARNING: Experiment found in proteinGroups,"
                          " but missing from experimental design:", exp_name)
                else:
                    if quantification_saint == "spc":
                        self.data[col] = self.data[col].astype(int)
                    else:
                        self.data[col] = self.data[col].astype(float)
                    self.quant_cols_saint.append(col)
                    if experimental_design.name2experiment[exp_name].attributes["Type"] == "T":
                        self.bait_cols.append(col)
            if col.startswith(self.quant_col_prefix_comppass):
                exp_name = col.split(self.quant_col_prefix_comppass)[-1]
                if exp_name in self.experimental_design.name2experiment:
                    if quantification_comppass == "spc":
                        self.data[col] = self.data[col].astype(int)
                    else:
                        self.data[col] = self.data[col].astype(float)
                    self.quant_cols_comppass.append(col)

        for name in self.experimental_design.name2experiment:
            if (self.quant_col_prefix_saint + name) not in self.quant_cols_saint or \
                    (self.quant_col_prefix_comppass + name) not in self.quant_cols_comppass:
                print("ERROR: Experiment found in experimental design, but missing from proteinGroups:", name)
                exit(1)

        print("Removed", sum(self.data["Reverse"] == "+"), "reverses")
        print("Removed", sum(self.data["Only identified by site"] == "+"), "only identified by site")
        print("Removed", sum(self.data["Potential contaminant"] == "+"), "potential contaminants")

        self.data = self.data[self.data["Reverse"] != "+"]
        self.data = self.data[self.data["Only identified by site"] != "+"]
        self.data = self.data[self.data["Potential contaminant"] != "+"]

        # remove proteins only found in controls
        print("Removed", sum(self.data[self.bait_cols].sum(axis=1) == 0), "proteins found only in controls")
        self.data = self.data[self.data[self.bait_cols].sum(axis=1) > 0]

        print("Kept", len(self.data.index), "proteins")

        self.data["Short protein IDs"] = self.data["Majority protein IDs"]
        self.data["Short Gene names"] = self.data["Gene names"]
        self.data.loc[self.data["Short Gene names"].isnull(), "Short Gene names"] = \
            self.data[self.data["Short Gene names"].isnull()]["Short protein IDs"]

        for index in self.data.index:
            split_prots = self.data.loc[index, "Short protein IDs"].split(";")
            split_genes = self.data.loc[index, "Short Gene names"].split(";")

            if len(split_prots) > 3:
                self.data.loc[index, "Short protein IDs"] = re.sub(" ", "-", ";".join(split_prots[0:3]) + ";+" + str(len(split_prots)-3))
            if len(split_genes) > 3:
                self.data.loc[index, "Short Gene names"] = re.sub(" ", "-", ";".join(split_genes[0:3]) + ";+" + str(len(split_genes)-3))

    def impute(self):
        1

    def correct_carry_over(self):
        1

    def write_prey_file(self, out_path):
        if self.quantification_saint == "spc":
            self.data[["Short protein IDs", "Sequence length", "Short Gene names"]].to_csv(out_path, index=False,
                                                                                           sep="\t", header=False)
        else:
            self.data[["Short protein IDs", "Short Gene names"]].to_csv(out_path, index=False, sep="\t", header=False)

    def write_bait_file(self, out_path):
        with open(out_path, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            for e in self.experimental_design.name2experiment.values():
                writer.writerow([e.attributes["Experiment Name"], e.attributes["Bait"], e.attributes["Type"]])

    def write_interaction_file(self, out_path):
        with open(out_path, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            for index, row in self.data.iterrows():
                for quant_col in self.quant_cols_saint:
                    exp_name = quant_col.split(self.quant_col_prefix_saint)[-1]
                    writer.writerow([exp_name,
                                    self.experimental_design.name2experiment[exp_name].attributes["Bait"],
                                    row["Short protein IDs"],
                                    row[quant_col]])

# CompPASS CSV File to be used in R
    def write_CompPASS(self, out_path):
        with open(out_path, 'w', newline='\n', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(["Experiment.ID", "Replicate", "Bait", "Prey", "Prey.Name", "Spectral.Count"])
            for index, row in self.data.iterrows():
                for quant_col in self.quant_cols_comppass:
                    exp_name = quant_col.split(self.quant_col_prefix_comppass)[-1]
                    exp_name_trim = re.match("(.+)_\d+$", exp_name).groups()[0]
                    if row[quant_col] > 0:
                        writer.writerow([self.experimental_design.name2experiment[exp_name].attributes["Bait"],
                                     self.experimental_design.name2experiment[exp_name].attributes["Replicate"],
                                     self.experimental_design.name2experiment[exp_name].attributes["Bait ID"],
                                     row["Short protein IDs"], row["Gene names"], row[quant_col]])

# Start ComPASS
    def to_CompPASS(self, out_path):
        self.write_CompPASS(os.path.join(out_path,"output/to_CompPASS.csv"))

    def to_SAINT(self, out_path):
        self.write_prey_file(os.path.join(out_path, "output", "prey.txt"))
        self.write_bait_file(os.path.join(out_path, "output/bait.txt"))
        self.write_interaction_file(os.path.join(out_path, "output/interaction.txt"))

    def to_SAINTq(self, out_path, normalize_control, compress_n_ctrl, compress_n_rep):
        self.write_toSAINTq_intensity_file(os.path.join(out_path, "intensities.txt"))
        self.write_toSAINTq_param_file("intensities.txt",
                                       os.path.join(out_path, "params.txt"),
                                       normalize_control, compress_n_ctrl, compress_n_rep)

    def write_toSAINTq_intensity_file(self, out_path):
        with open(out_path, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')

            for index, row in self.data.iterrows():
                if index == 0:
                    header1 = [""]
                    header2 = [""]
                    header3 = ["Protein"]
                    for quant_col in self.quant_cols_saint:
                        exp_name = quant_col.split(self.quant_col_prefix_saint)[-1]
                        header1.append(self.experimental_design.name2experiment[exp_name].attributes["Type"])
                        header2.append(self.experimental_design.name2experiment[exp_name].attributes["Bait"])
                        header3.append(exp_name)

                    writer.writerow(header1)
                    writer.writerow(header2)
                    writer.writerow(header3)

                quant_row = [row["Short protein IDs"]]
                for quant_col in self.quant_cols_saint:
                    quant_row.append(row[quant_col])
                writer.writerow(quant_row)

    def write_toSAINTq_param_file(self, intensities_path, out_path, normalize_control, compress_n_ctrl, compress_n_rep):
        outfile = open(out_path, 'w')

        outfile.write("input_filename=" + intensities_path + "\n")

        outfile.write("input_level=protein\n")

        outfile.write("protein_colname=Protein\n")

        outfile.write("normalize_control=")
        if normalize_control==1:
            outfile.write("true\n")
        else:
            outfile.write("false\n")

        outfile.write("compress_n_ctrl=")
        outfile.write(str(compress_n_ctrl) + "\n")

        outfile.write("compress_n_rep=")
        outfile.write(str(compress_n_rep) + "\n")
        outfile.close()

    def calc_CompPASS(self):
        # num baits
        bait2protID = dict()
        for e in self.experimental_design.name2experiment.values():
            bait2protID[e.attributes["Bait"]] = e.attributes["Bait ID"]
        k = len(bait2protID)

        # for each prey, bait: compute total quant, replicate count, average quant
        prey2bait2tot = defaultdict(dict)
        prey2bait2rep = defaultdict(dict)
        prey2avg = defaultdict(float)
        for index, row in self.data.iterrows():
            prey = row["Short protein IDs"]
            for e in self.experimental_design.name2experiment.values():
                bait = e.attributes["Bait"]
                exp_name = e.attributes["Experiment Name"]

                if bait not in prey2bait2tot[prey]:
                    prey2bait2tot[prey][bait] = 0
                    prey2bait2rep[prey][bait] = 0
                if self.quantification_comppass == "spc":
                    prey2bait2tot[prey][bait] += row[self.quant_col_prefix_comppass + exp_name]
                else:
                    prey2bait2tot[prey][bait] += log2(1+row[self.quant_col_prefix_comppass + exp_name])

                if row[self.quant_col_prefix_comppass + exp_name] > 0:
                    prey2bait2rep[prey][bait] += 1

                if self.quantification_comppass == "spc":
                    prey2avg[prey] += row[self.quant_col_prefix_comppass + exp_name] / len(self.experimental_design.name2experiment)
                else:
                    prey2avg[prey] += log2(1+row[self.quant_col_prefix_comppass + exp_name]) / len(self.experimental_design.name2experiment)

        # for each prey: compute frequency
        prey2freq = dict()
        for prey in prey2bait2tot:
            tot_bait = 0
            for bait in prey2bait2tot[prey]:
                if prey2bait2tot[prey][bait] > 0 and prey != bait2protID[bait]:
                    tot_bait += 1
            prey2freq[prey] = tot_bait

        # for each prey: compute std
        prey2stdev = defaultdict(float)
        prey2W = dict()
        for index, row in self.data.iterrows():
            prey = row["Short protein IDs"]
            for e in self.experimental_design.name2experiment.values():
                exp_name = e.attributes["Experiment Name"]
                if self.quantification_comppass == "spc":
                    prey2stdev[prey] += pow(row[self.quant_col_prefix_comppass + exp_name] - prey2avg[prey], 2)
                else:
                    prey2stdev[prey] += pow(log2(1 + row[self.quant_col_prefix_comppass + exp_name]) - prey2avg[prey], 2)

                prey2stdev[prey] = sqrt(prey2stdev[prey] / len(self.experimental_design.name2experiment))

            if prey2avg[prey] > 0:
                prey2W[prey] = max(1.0, prey2stdev[prey] / prey2avg[prey])
            else:
                prey2W[prey] = 0.0

        # calculate WD-score
        prey2bait2scores = defaultdict(dict)
        for index, row in self.data.iterrows():
            prey = row["Short protein IDs"]

            for bait in bait2protID:
                if prey2avg[prey] > 0:
                    wd = sqrt(pow(((k / max(1,prey2freq[prey])) * prey2W[prey]), prey2bait2rep[prey][bait]) * prey2bait2tot[prey][bait])
                    d = sqrt(pow((k / max(1,prey2freq[prey])), prey2bait2rep[prey][bait]) * prey2bait2tot[prey][bait])
                else:
                    wd = 0.0
                    d = 0.0

                prey2bait2scores[prey][bait] = [str(wd), str(d), str(prey2freq[prey])]

        return prey2bait2scores

#    def align_scores(self, saint_path, prey2bait2comppass, out_path):
#        with open(saint_path, encoding='utf-8-sig') as csvfile:
#            reader = csv.reader(csvfile, delimiter="\t")
#            header = next(reader)

#            with open(out_path, 'w') as csvfile_out:
#                writer = csv.writer(csvfile_out, delimiter='\t')

#                new_header = header
#                new_header.extend(["WD-score", "D-score", "frequency"])
#                writer.writerow(new_header)
#                for row in reader:
#                    prey = row[1]
#                    bait = row[0]
#                    new_row = row

#                    new_row.extend(prey2bait2comppass[prey][bait])
#                    writer.writerow(new_row)





