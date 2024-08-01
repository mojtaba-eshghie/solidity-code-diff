import pickle
import pprint
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def load_pickles():
    difft_file = open("results-difft.pickle", "rb")
    gt_file = open("results-GT.pickle", "rb")

    difft = pickle.load(difft_file)
    GT = pickle.load(gt_file)

    return {"GT": GT, "difft": difft}


def setup(n_mut):
    res = [0] * n_mut
    count = [0] * n_mut
    mut_res = {}
    mut_count = {}

    return {"res": res, "count": count, "mut_res": mut_res, "mut_count": mut_count}


def remove_mut_operators(diff_results, operators):
    for diff_tool in diff_results:
        for contract in diff_results[diff_tool]:
            for i in range(len(diff_results[diff_tool][contract])):
                for r in operators:
                    if r in diff_results[diff_tool][contract][i].keys():
                        diff_results[diff_tool][contract][i].pop(r)


def analyze_diffs(diffs, res_dict, n_mut):
    is_GT = isinstance(next(iter(next(iter(diffs.values()))[0].values())), list)

    for contract in diffs:
        for i in range(len(diffs[contract])):
            for mut in diffs[contract][i]:
                if not diffs[contract][i][mut] == []:
                    res_dict["count"][i] += 1
                    if is_GT:
                        res_dict["res"][i] += diffs[contract][i][mut][0]
                    else:
                        res_dict["res"][i] += diffs[contract][i][mut]
                    
                    if not mut in res_dict["mut_res"].keys():
                        res_dict["mut_res"][mut] = [0] * n_mut
                        res_dict["mut_count"][mut] = [0] * n_mut 

                    if is_GT:
                        res_dict["mut_res"][mut][i] += diffs[contract][i][mut][0]
                    else:
                        res_dict["mut_res"][mut][i] += diffs[contract][i][mut] 
                    res_dict["mut_count"][mut][i] += 1    

    #Calculate average results 
    res_dict["res"] = [i / j for i, j in zip( res_dict["res"],  res_dict["count"])]

    for mut in res_dict["mut_res"]:
        for i in range(len(np.trim_zeros( res_dict["mut_res"][mut]))):
            res_dict["mut_res"][mut] = np.trim_zeros( res_dict["mut_res"][mut])
            res_dict["mut_count"][mut] = np.trim_zeros( res_dict["mut_count"][mut])
            res_dict["mut_res"][mut][i] /= res_dict["mut_count"][mut][i]


def print_summary(GT, difft):
    print("====================================")
    print("difftastic average: " + str(difft["res"]))
    print("successful # of mutations: ", difft["count"])
    print("Gumtree average: " + str(GT["res"]))
    print("successful # of mutations: ", GT["count"])
    print("====================================")


    plt.bar(x = [1,2,3,4,5,6,7,8,9,10], height = difft["count"], color = 'b')
    plt.xticks([1,2,3,4,5,6,7,8,9,10])
    plt.title("Number of Times X Amount of Mutations Could be Applied to Files")
    plt.xlabel("Number of Mutations")
    plt.ylabel("Successful Mutations")
    plt.show()


def print_by_mutation(res_dict):
    
    print("EDIT SCRIPT LENGTH BY MUTATION:")
    print("=====================================================")
    mut_sorted = sorted(res_dict["mut_res"].items(), key = lambda x: x[1][len(x[1])-1])
    for mut in mut_sorted:
        print(mut[0], mut[1], " counts:", res_dict["mut_res"][mut[0]])
    print("=====================================================")


def save_as_csv(filename):
    with open(filename, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter= " ", quotechar='|')
        for key in GT_res_dict["mut_res"]:
            writer.writerow([key] + GT_res_dict["mut_res"][key])
        writer.writerow("           ")
        for key in difft_res_dict["mut_res"]:
            writer.writerow([key] + difft_res_dict["mut_res"][key])


def scatter_with_avg_plot(GT, difft):
    plt.plot(range(10) ,difft["res"], label=("difft_avg"), color="blue")
    for mut in difft["mut_res"]:
        x = [i+0.985 for i in range(len(difft["mut_res"][mut]))]
        plt.scatter(x, difft["mut_res"][mut], color="blue", s=6)

    plt.plot(range(10), GT["res"], label ="GT_avg", color = "red")
    for mut in GT["mut_res"]:
        x = [i+1.015 for i in range(len(GT["mut_res"][mut]))]
        plt.scatter(x, GT["mut_res"][mut], color="red", s=6)

    plt.title("Average Edit Distance per Mutation and Tool")
    plt.ylabel("edit actions")
    plt.xlabel("# of mutations")
    plt.minorticks_on()
    plt.xticks(range(11))
    plt.legend()
    plt.show()


def box_plot(data, offset):
    box(data[0], offset, "green")
    box(data[1], -offset, "blue")

    green_patch = mpatches.Patch(color="green", label="Gumtree")
    blue_patch = mpatches.Patch(color="blue", label="difftastic")
    plt.legend(handles=[green_patch, blue_patch])

    x = [1,2,3,4,5,6,7,8,9,10]
    plt.ylabel("edit actions")
    plt.xlabel("# of mutations")
    plt.xticks(x, labels=x)
    plt.title("Edit Distances")
    plt.show()


def box(data, offset, color):
    d = [[] for _ in range(10)]
    for key in data.keys():
        for i in range(len(data[key])):
            d[i].append(data[key][i])
    
    x = [x - offset for x in [1,2,3,4,5,6,7,8,9,10]]
   
    plt.boxplot(d, widths=0.27, patch_artist=True,
        positions=x,
        showmeans=True, showfliers=False,
        medianprops={"color": "white", "linewidth": 1},
        boxprops={"facecolor":  color, "edgecolor": "black",
                  "linewidth": 0.25},
        whiskerprops={"color": color, "linewidth": 1},
        meanprops={"markerfacecolor": color, "markeredgecolor": "black"},
        capprops={"color": color, "linewidth": 1})


def bar_by_mut_plot(data, offset, oper, title):
    
    fig, ax = plt.subplots(layout='constrained')

    bar(data[0], -offset, "green", ax, oper)
    bar(data[1], offset, "blue", ax, oper)

    red_patch = mpatches.Patch(color="lime", label="Gumtree")
    blue_patch = mpatches.Patch(color="b", label="difftastic")
    plt.legend(handles=[red_patch, blue_patch])
    
    plt.title(title)
    plt.show()


def bar(data, offset, color, subplot, oper):  
    #greens = ["darkgreen", "green", "forestgreen", "seagreen", "lime", "lawngreen", "lightgreen", "chartreuse", "greenyellow", "springgreen"]
    greens = ['#173d00','#225009','#2e6413','#3a781e','#478d29','#53a334','#60b940','#6cd04c','#79e759','#85ff66']
    #blues= ["navy", "midnightblue", "darkblue", "mediumblue", "blue","royalblue",  "deepskyblue", "turquoise","cyan", "paleturquoise"]
    blues = ['#140052','#071f70','#003b8b','#0056a2','#0071b6','#008dca','#00aada','#0bc6e8','#4ae3f4','#75ffff']

    colors = "none"
    if color == "blue":
        colors = blues
    elif color == "green":
        colors = greens

    x = np.arange(len(oper))
    n = 0
    for key in oper:
        plt_data = np.around(data[key], 2)
        for i in range(len(plt_data)-1, -1, -1):
            subplot.bar(x = x[n] + offset, width =  0.25, height = plt_data[i], color = colors[i])
        n += 1
    subplot.set_xticks(x, oper)
    

def calc_corr(data):
    x = [1,2,3,4,5,6,7,8,9,10]
    count = [0 for _ in range(10)]
    
    d = [0 for _ in range(10)]
    for key in data.keys():
        for i in range(len(data[key])):
            d[i] += data[key][i]
            count[i] += 1
    
    res = [i / j for i, j in zip(d, count)]

    r = np.corrcoef(x, res)
    r2 = np.polyfit(x, res, 1)
    print("Slope: ", r2, "\nCorrelation:", r, '\n')        


if __name__ ==  '__main__':
    num_mut = 10
    pickles = load_pickles()

    remove_mut_operators(pickles, ["AVR","SCEC"])

    GT_res_dict = setup(num_mut)
    analyze_diffs(pickles["GT"], GT_res_dict, num_mut)

    difft_res_dict = setup(num_mut)
    analyze_diffs(pickles["difft"], difft_res_dict, num_mut)

    #print_summary(GT_res_dict, difft_res_dict)
    #print_by_mutation(difft_res_dict)
    #print_by_mutation(GT_res_dict)

    #scatter_with_avg_plot(GT_res_dict, difft_res_dict)
    box_plot((GT_res_dict["mut_res"], difft_res_dict["mut_res"]), 0.15)
   

    x = [1,2,3,4,5,6,7,8,9,10]
    r = np.corrcoef(x, GT_res_dict["res"])
    r2 = np.polyfit(x, GT_res_dict["res"], 1)
    print("Gumtree slope: ", r2, "\nGumtree correlation:", r, '\n')
    calc_corr( GT_res_dict["mut_res"])

    
    r = np.corrcoef(x, difft_res_dict["res"])
    r2 = np.polyfit(x, difft_res_dict["res"], 1)
    print("difft slope: ", r2, "\ndifft correlation:", r)
    calc_corr( difft_res_dict["mut_res"])
   
    

    # removed operators due to errors in testing: "AVR",  "SCEC"
    opers = []
    opers.append( ["BLR", "HLR", "ILR", "SLR"])
    opers.append( ["AOR", "BOR", "DOD", "ECS", "ICM", "MCR", "UORD","VUR"])
    opers.append( ["CBD", "CCD", "CSC", "EED", "EHC", "OLFD", "ORFD", "RSD"])
    opers.append( ["ACM", "LSC", "MOC", "MOD", "MOI", "MOR", "RVS"])
    opers.append( ["BCRD","DLR","ER","ETR","FVR","GVR","PKD","SFR","SKD","SKI","TOR","VVR"])
    titles = ["Mutated Literals", "Mutated Operators & Type Specifications", "Mutated Code Blocks", "Mutated Arguments & Modifers", "Other Mutations"]

    for i in range(len(opers)):
        bar_by_mut_plot((GT_res_dict["mut_res"], difft_res_dict["mut_res"]), 0.15, opers[i], titles[i])
    

    
    