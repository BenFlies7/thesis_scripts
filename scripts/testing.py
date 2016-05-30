import numpy as np
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import pandas as pd

# Read data
data = pd.read_excel(input_file, sheetname=sheet)

# Create three sets of the lists to be compared
set_1 = set(data[compare[0]].dropna())
set_2 = set(data[compare[1]].dropna())
set_3 = set(data[compare[2]].dropna())

# Create a third set with all elements of the two lists
union = set_1.union(set_2).union(set_3)

# Gather names of all elements and list them in groups
lists = [[], [], [], [], [], [], []]
for gene in union:
    if (gene in set_1) and (gene not in set_2) and (gene not in set_3):
        lists[0].append(gene)
    elif (gene in set_1) and (gene in set_2) and (gene not in set_3):
        lists[1].append(gene)
    elif (gene in set_1) and (gene not in set_2) and (gene in set_3):
        lists[2].append(gene)
    elif (gene in set_1) and (gene in set_2) and (gene in set_3):
        lists[3].append(gene)
    elif (gene not in set_1) and (gene in set_2) and (gene not in set_3):
        lists[4].append(gene)
    elif (gene not in set_1) and (gene in set_2) and (gene in set_3):
        lists[5].append(gene)
    elif (gene not in set_1) and (gene not in set_2) and (gene in set_3):
        lists[6].append(gene)

# Write gene lists to file
ew = pd.ExcelWriter('../Gene lists/Venn lists/' + compare[0] + ' & '
                    + compare[1] + ' & ' + compare[2] + ' gene lists.xlsx')

pd.DataFrame(lists[0], columns=[compare[0]]) \
    .to_excel(ew, sheet_name=compare[0], index=False)

pd.DataFrame(lists[1], columns=[compare[0] + ' & ' + compare[1]]) \
    .to_excel(ew, sheet_name=compare[0] + ' & ' + compare[1], index=False)

pd.DataFrame(lists[2], columns=[compare[0] + ' & ' + compare[2]]) \
    .to_excel(ew, sheet_name=compare[0] + ' & ' + compare[2], index=False)

pd.DataFrame(lists[3], columns=['All']) \
    .to_excel(ew, sheet_name='All', index=False)

pd.DataFrame(lists[4], columns=[compare[1]]) \
    .to_excel(ew, sheet_name=compare[1], index=False)

pd.DataFrame(lists[5], columns=[compare[1] + ' & ' + compare[2]]) \
    .to_excel(ew, sheet_name=compare[1] + ' & ' + compare[2], index=False)

pd.DataFrame(lists[6], columns=[compare[2]]) \
    .to_excel(ew, sheet_name=compare[2], index=False)

ew.save()

# Count the elements in each group
subsets = [len(lists[0]), len(lists[4]), len(lists[1]), len(lists[6]),
           len(lists[2]), len(lists[5]), len(lists[3])]

# Basic venn diagram
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1)
v = venn3(subsets, (compare[0], compare[1], compare[2]), ax=ax)
c = venn3_circles(subsets)

# Annotation
ax.annotate('Total genes:\n' + str(len(union)),
            xy=v.get_label_by_id('111').get_position() - np.array([-0.5,
                                                                   0.05]),
            xytext=(0,-70), ha='center', textcoords='offset points',
            bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.3))

# Title
plt.title(compare[0] + ' & ' + compare[1] + ' & ' + compare[2] +
          ' gene expression overlap')

plt.show()
