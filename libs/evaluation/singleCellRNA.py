import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib.colors as colors
import numpy as np
import pandas as pd

#used

def single_cell_RNA(singleCell, RNAseq, attr):
    attribute_values, fold_changes, color = [], [], []
    name = []
    for seq in singleCell.sequences:
        if seq.gene in list(RNAseq.sequences_by_name.keys()):
            attribute_values.append(getattr(seq, attr))
            fold_changes.append(RNAseq.sequences_by_name[seq.gene].fold_change)
            color.append(RNAseq.sequences_by_name[seq.gene].p_value)
            name.append(seq.gene)

    sc = plt.scatter(attribute_values, fold_changes, c=color, norm=colors.LogNorm(), s=20, cmap="coolwarm")
    plt.colorbar(sc).set_label("log p value")

    m, b, r_value, p_value, std_err = linregress(attribute_values, fold_changes)
    x = np.linspace(np.amin(attribute_values), np.amax(attribute_values), 10)
    plt.plot(x, m * x + b, color="black",
             label="slope: " + str(np.round(m, decimals=3)) + ", intersection:" + str(np.round(b, decimals=1))
                   + ", r_value: " + str(np.round(r_value, decimals=4)) + ", p_values: " + str(p_value))
    plt.legend(bbox_to_anchor=(1.3, 1.2))

    plt.xlabel(singleCell.name+"_"+attr)
    plt.ylabel(RNAseq.name)
    plt.subplots_adjust(top=0.8)
    plt.savefig("results/images/singleCellvsRNA___" + attr + "___" + singleCell.name + "___" + RNAseq.name + ".png")
    plt.close()
    pd.DataFrame(np.array([name, attribute_values, fold_changes, color]).T).to_csv(
        "results/data/singleCellvsRNA___" + attr + "___" + singleCell.name + "___" + RNAseq.name + ".csv", index=False, header=False)