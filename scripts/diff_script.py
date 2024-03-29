#Runs gumtree to get data on diffs, puts result in 2d matrix with axes contract and number of mutations.
#Stores results in pickle file. 

# -----------Gumtree Results----------------
# Contracts |      Diff Results/# Mutations
# ------------------------------------------
# Contract1 | 1 mutation, 2 mutations, ...
# Contract2 | 1 mutation, ...
# ...
# Contractn | ...
# ------------------------------------------

import sys
import os
import xml.etree.ElementTree as ET
import subprocess
import pickle
import json
import time
import concurrent.futures as cc


def parse_input():
    if len(sys.argv) != 4:
        raise Exception("Please provide arguments in the form: [PATH TO MUTANTS], [DIFFING TOOL] (GT/difft), [LANGUAGE] (Solidity/Python)!")
    
    contracts_path = sys.argv[1]
    diff_tool = sys.argv[2]
    file_ending = ""
    if sys.argv[3] == "Solidity":
        file_ending = ".sol"
    elif sys.argv[3] == "Python":
        file_ending = ".py"
    else:
        raise Exception("Invalid programming language")

    return contracts_path, diff_tool, file_ending

#Uses gumtree to get the diff data of two files
def get_GT_diff_data(filepath1, filepath2):
    save_full_diff = False

    diff = subprocess.check_output('gumtree textdiff -f XML ' +  filepath1 + " " + filepath2, shell=True).decode()
    if not diff:
        print("mutant causing error:" + filepath2)
        return False
    
    #wrap result to get single XML root and convert to tree
    diff = diff.split('\n', 1)
    diff = diff[0] + "<X>" + diff[1] + "</X>"
    tree = ET.fromstring(diff)

    #Get number of edit actions and matches from tree
    n_edits = len(tree.findall('actions')[0])
    n_matches = len(tree.findall('matches')[0])
    
    #Append full diff to res if flag is set and return res
    res = [n_edits, n_matches]
    if save_full_diff:
        res.append(diff)
    return res

#Uses difftastic to get the diff data
def get_diffts_data(filepath1, filepath2):
    os.environ['DFT_UNSTABLE'] = 'yes'
    diff = subprocess.check_output('difft --display json ' +  filepath1 + " " + filepath2, shell=True).decode()
    if not diff:
        return False 
    
    diff = json.loads(diff)
    if diff["status"] == "unchanged":
        print("WARNING: difft failed to detect change!!!")
        return False
    
    count = 0
    for a in diff["chunks"]:
        count += len(a)
    return count
    
#Gets the diffs between all mutants of a contract in a directory
def get_contract_diffs(contract_path, contract, file_ending, result, diff_tool):
    res = []
    unmutated_path = contract_path + "/" + "original/" + contract + file_ending
    num_contracts = len(os.listdir(contract_path))
    for i in range(1, num_contracts):
        operators = os.listdir(contract_path + "/" + str(i))
        diffs = {}
        for op in operators:
            mutated_path = contract_path + "/" + str(i) + "/" + op + "/" + contract + file_ending
            if diff_tool == "GT":
                diff = get_GT_diff_data(unmutated_path, mutated_path)
                if not diff:
                    diff = []
            elif diff_tool == "difft":
                diff = get_diffts_data(unmutated_path, mutated_path)
                if not diff:
                    diff = 0
            else:
                raise Exception("Error: invalid diff tool provided!")
            diffs[op] = diff
        res.append(diffs)
    result[contract] = res

        
#Returns complete 2d matrix containing diff data for all mutants. Calculates each contract's diff concurrently.
def calculate_diffs(contracts_path, file_ending, diff_tool):
    contracts = os.listdir(contracts_path)
    res = {}
    executor = cc.ThreadPoolExecutor(max_workers = os.cpu_count())
    futures = []

    for c in contracts:
        futures.append(executor.submit(get_contract_diffs, contracts_path+c, c, file_ending, res, diff_tool))

    while len(futures) > 0:
        time.sleep(1)
        for i in range(len(futures)-1, -1, -1):
            if(futures[i].done()):
                futures.pop(i)
        print('Contracts done: ' + str(len(contracts) - len(futures)) + '/' + str(len(contracts)), end=('\r'))
    return res

#Saves results as a python object in a results file
def save_res_to_file(results, diff_tool):
    out = "./results/results" + "-" + diff_tool +  ".pickle"
    print("\nDumping results to " + out)
    out_file = open(out, "wb")
    pickle.dump(results, out_file)

if __name__ ==  '__main__':
    seconds = int(time.time())
    
    contracts_path, diff_tool, file_ending = parse_input()
    res = calculate_diffs(contracts_path, file_ending, diff_tool)
    save_res_to_file(res, diff_tool)
    
    m, s = divmod(int(time.time()) -  seconds, 60)
    print("Generated diffs in " + str(m) + "m:" + str(s) + "s")


    
    