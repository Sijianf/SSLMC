
import os
import numpy as np
from collections import defaultdict

''' 
load_data_from_file, get_drugs_targets_names: 
    they are modified (ugly) by sijian in 2024.10.30 to be compatible for simulation study 
'''

def load_data_from_file(dataset, folder):
    if dataset == "simulation":
        # Read in the "Y.txt" file without skipping the first line
        with open(os.path.join(folder, "Y.txt"), "r") as inf:
            int_array = [line.strip("\n").split() for line in inf]
        
        # Convert int_array to a NumPy array for the interaction matrix
        intMat = np.array(int_array, dtype=np.float64)
        
        # Create drug_sim and target_sim as identity matrices based on Y's dimensions
        drug_sim = np.identity(intMat.shape[0], dtype=np.float64)
        target_sim = np.identity(intMat.shape[1], dtype=np.float64)
    else:
        # For non-simulation datasets
        with open(os.path.join(folder, dataset+"_admat_dgc.txt"), "r") as inf:
            # next(inf)  # skip the first line if needed
            int_array = [line.strip("\n").split() for line in inf]

        with open(os.path.join(folder, dataset+"_simmat_dc.txt"), "r") as inf:
            # next(inf)  # skip the first line if needed
            drug_sim = [line.strip("\n").split() for line in inf]

        with open(os.path.join(folder, dataset+"_simmat_dg.txt"), "r") as inf:
            # next(inf)  # skip the first line if needed
            target_sim = [line.strip("\n").split() for line in inf]

        # Convert lists to NumPy arrays for the matrices
        intMat = np.array(int_array, dtype=np.float64)      # drug-target interaction matrix
        drug_sim = np.array(drug_sim, dtype=np.float64)     # drug similarity matrix
        target_sim = np.array(target_sim, dtype=np.float64) # target similarity matrix

    return intMat, drug_sim, target_sim


def get_drugs_targets_names(dataset, folder):
    if dataset == "simulation":
        # Read the "Y.txt" file without skipping the first line
        with open(os.path.join(folder, "Y.txt"), "r") as inf:
            int_array = [line.strip("\n").split() for line in inf]
        
        # Assign drug and target names based on the dimensions of int_array
        num_drugs = len(int_array)         # Number of rows in Y
        num_targets = len(int_array[0])    # Number of columns in Y
        
        # Assign drugs as "d1", "d2", ..., "d<num_drugs>"
        drugs = [f"d{i+1}" for i in range(num_drugs)]
        
        # Assign targets as "t1", "t2", ..., "t<num_targets>"
        targets = [f"t{j+1}" for j in range(num_targets)]
    
    else:
        # For non-simulation datasets
        with open(os.path.join(folder, dataset+"_admat_dgc.txt"), "r") as inf:
            int_array = [line.strip("\n").split() for line in inf]
        
        # Assign drug and target names based on the dimensions of int_array
        num_drugs = len(int_array)         # Number of rows in data
        num_targets = len(int_array[0])    # Number of columns in data
        
        # Assign drugs as "d1", "d2", ..., "d<num_drugs>"
        drugs = [f"d{i+1}" for i in range(num_drugs)]
        
        # Assign targets as "t1", "t2", ..., "t<num_targets>"
        targets = [f"t{j+1}" for j in range(num_targets)]    
            
        # Use following if your file have row and column names
        # with open(os.path.join(folder, dataset+"_admat_dgc.txt"), "r") as inf:
        #     drugs = [line.strip("\n").split()[0] for line in inf]
        #     targets = next(inf).strip("\n").split() 

    return drugs, targets


def cross_validation(intMat, seeds, cv=0, num=10):
    cv_data = defaultdict(list)
    for seed in seeds:
        num_drugs, num_targets = intMat.shape
        prng = np.random.RandomState(seed)
        if cv == 0:
            index = prng.permutation(num_drugs)
        if cv == 1:
            index = prng.permutation(intMat.size)
        # Use integer division
        step = index.size // num
        for i in range(num):
            if i < num-1:
                ii = index[i*step:(i+1)*step]
            else:
                ii = index[i*step:]
            if cv == 0:
                test_data = np.array([[k, j] for k in ii for j in range(num_targets)], dtype=np.int32)
            elif cv == 1:
                test_data = np.array([[k/num_targets, k % num_targets] for k in ii], dtype=np.int32)
            x, y = test_data[:, 0], test_data[:, 1]
            test_label = intMat[x, y]
            W = np.ones(intMat.shape)
            W[x, y] = 0
            cv_data[seed].append((W, test_data, test_label))
    return cv_data


def train(model, cv_data, intMat, drugMat, targetMat):
    aupr, auc = [], []
    for seed in cv_data.keys():
        for W, test_data, test_label in cv_data[seed]:
            model.fix_model(W, intMat, drugMat, targetMat, seed)
            aupr_val, auc_val = model.evaluation(test_data, test_label)
            aupr.append(aupr_val)
            auc.append(auc_val)
    return np.array(aupr, dtype=np.float64), np.array(auc, dtype=np.float64)


def svd_init(M, num_factors):
    from scipy.linalg import svd
    U, s, V = svd(M, full_matrices=False)
    ii = np.argsort(s)[::-1][:num_factors]
    s1 = np.sqrt(np.diag(s[ii]))
    U0, V0 = U[:, ii].dot(s1), s1.dot(V[ii, :])
    return U0, V0.T


def mean_confidence_interval(data, confidence=0.95):
    import scipy as sp
    import scipy.stats
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, h


def write_metric_vector_to_file(auc_vec, file_name):
    np.savetxt(file_name, auc_vec, fmt='%.6f')


def load_metric_vector(file_name):
    return np.loadtxt(file_name, dtype=np.float64)
