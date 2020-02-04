import pandas as pd
import os

class ProteinGroups:
    def __init__(self, experimental_design, file_path, quantification):
        print("\nParsing MaxQuant proteinGroups: " + file_path)

        self.quantification = quantification
        self.experimental_design = experimental_design
        self.data = pd.read_csv(file_path, sep="\t", quotechar="'")

        if quantification == "LFQ":
            quant_col_prefix = "LFQ intensity "
        elif quantification == "intensity":
            quant_col_prefix = "Intensity "
        elif quantification == "spc":
            quant_col_prefix = "MS/MS count "

        self.quant_cols = []
        self.bait_cols = []

        # make sure experimental design matches the data
        for col in self.data.columns:
            if col.startswith(quant_col_prefix):
                exp_name = col.split(quant_col_prefix)[-1]
                if exp_name not in self.experimental_design.name2experiment:
                    print("WARNING: Experiment found in proteinGroups,"
                          " but missing from experimental design:", exp_name)
                else:
                    if quantification == "spc":
                        self.data[col] = self.data[col].astype(int)
                    else:
                        self.data[col] = self.data[col].astype(float)
                    self.quant_cols.append(col)
                    if experimental_design.name2experiment[exp_name].attributes["Type"] == "T":
                        self.bait_cols.append(col)

        for name in self.experimental_design.name2experiment:
            if (quant_col_prefix + name) not in self.quant_cols:
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
                self.data.loc[index, "Short protein IDs"] = ";".join(split_prots[0:3]) + ";+" + str(len(split_prots)-3)
            if len(split_genes) > 3:
                self.data.loc[index, "Short Gene names"] = ";".join(split_genes[0:3]) + ";+" + str(len(split_genes)-3)

    def impute(self):
        1

    def correct_carry_over(self):
        1

    def write_prey_file(self, out_path):
        self.data[["Short protein IDs", "Sequence length", "Short Gene names"]].to_csv(out_path, index=False, sep="\t", header=False)

    def write_bait_file(self, out_path):
        1

    def write_interaction_file(self, out_path):
        1

    def to_SAINTv2(self, out_path):
        self.write_prey_file(os.path.join(out_path, "prey.txt"))
        self.write_bait_file(out_path)
        self.write_interaction_file(out_path)

    def to_SAINTexpress(self, out_path):
        1


