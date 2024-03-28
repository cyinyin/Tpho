"""
Tpho: A Python toolkit  to predicting the optimal catalytic temperature of phosphatase
"""

# Imports
# ============#
import pandas as pd
import os
from collections import Counter
from tqdm import tqdm
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np

# import PseudoAAC
# from PyPro import GetProDes
# from CTD import CalculateD
#
# from unit.path import Args
from sklearn.neighbors import KNeighborsRegressor
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from predicate_phosphatase_temperature import resreg
import sklearn.metrics as sm
from sklearn.feature_selection import SelectKBest, f_regression

# Variables
# =================#

aalist = list('ACDEFGHIKLMNPQRSTVWY')
AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

_Hydrophobicity = {"A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85, "E": -0.74, "G": 0.48,
                   "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12, "S": -0.18,
                   "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08}

_hydrophilicity = {"A": -0.5, "R": 3.0, "N": 0.2, "D": 3.0, "C": -1.0, "Q": 0.2, "E": 3.0, "G": 0.0, "H": -0.5,
                   "I": -1.8, "L": -1.8, "K": 3.0, "M": -1.3, "F": -2.5, "P": 0.0, "S": 0.3, "T": -0.4, "W": -3.4,
                   "Y": -2.3, "V": -1.5}

_residuemass = {"A": 15.0, "R": 101.0, "N": 58.0, "D": 59.0, "C": 47.0, "Q": 72.0, "E": 73.0, "G": 1.000, "H": 82.0,
                "I": 57.0, "L": 57.0, "K": 73.0, "M": 75.0, "F": 91.0, "P": 42.0, "S": 31.0, "T": 45.0, "W": 130.0,
                "Y": 107.0, "V": 43.0}

_pK1 = {"A": 2.35, "C": 1.71, "D": 1.88, "E": 2.19, "F": 2.58, "G": 2.34, "H": 1.78, "I": 2.32, "K": 2.20, "L": 2.36,
        "M": 2.28, "N": 2.18, "P": 1.99, "Q": 2.17, "R": 2.18, "S": 2.21, "T": 2.15, "V": 2.29, "W": 2.38, "Y": 2.20}

_pK2 = {"A": 9.87, "C": 10.78, "D": 9.60, "E": 9.67, "F": 9.24, "G": 9.60, "H": 8.97, "I": 9.76, "K": 8.90, "L": 9.60,
        "M": 9.21, "N": 9.09, "P": 10.6, "Q": 9.13, "R": 9.09, "S": 9.15, "T": 9.12, "V": 9.74, "W": 9.39, "Y": 9.11}

_pI = {"A": 6.11, "C": 5.02, "D": 2.98, "E": 3.08, "F": 5.91, "G": 6.06, "H": 7.64, "I": 6.04, "K": 9.47, "L": 6.04,
       "M": 5.74, "N": 10.76, "P": 6.30, "Q": 5.65, "R": 10.76, "S": 5.68, "T": 5.60, "V": 6.02, "W": 5.88, "Y": 5.63}
this_dir, this_filename = os.path.split(__file__)
singledatalist = []
alldatalist = []
namelist= []
seqlist = []
toptlist = []
labelslist = []
instability_indexs = []
molecular_weights = []
aromaticities = []
charge_at_pHs = []
flexibilities = []
isoelectric_points = []
secondaries = []
structures = []
fractions = []
gravies = []
protein_scales = []
selected_indexs = []
X_feature = []
bins = [30, 50, 65, 85]
sel = SelectKBest(score_func=f_regression, k=50)


# Get pptopt model
# ======================#


def do_count(seq):
    """
        Count the number of dipeptide in the sequence
    """
    dimers = Counter()
    for i in range(len(seq) - 1):
        dimers[seq[i:i + 2]] += 1.0
    return dimers


def getAAC(seq):
    """
        Calculate amino acid composition for a protein sequences

        Parameters
        -----------
        seq: str
        Sequence of phosphatase

        Returns
        ----------
        aac_array
        A array of amino acid frequency
        """

    aac = []
    for x in aalist:
        aac.append(seq.count(x) / len(seq))
    return np.array(aac)


def get_dimer_frequency(seq):
    """
        Calculating the dipeptide frequency of proteins

        Parameters
        -----------
        seq: str
        Sequence of phosphatase

        Returns
        ----------
        aac_array
        A array of dipeptide frequency of proteins
    """
    amino_acids = ['M', 'I', 'K', 'V', 'L', 'D', 'G', 'T', 'S', 'F', 'E', 'H', 'Q', 'A', 'R', 'P', 'Y',
                   'N', 'C', 'W']
    dimers_fq = dict()
    seq_matrix = np.zeros([20, 20])
    for a1 in amino_acids:
        for a2 in amino_acids:
            dimers_fq[a1 + a2] = do_count(seq).get(a1 + a2, 0.0)
    number_of_aa_in_fasta = sum(dimers_fq.values())
    if number_of_aa_in_fasta != 0:
        for key, value in dimers_fq.items():
            dimers_fq[key] = value / number_of_aa_in_fasta
    for index, item in enumerate(amino_acids):
        for inde, ite in enumerate(amino_acids):
            seq_matrix[index][inde] = dimers_fq[item + ite]
    return np.array(np.squeeze(seq_matrix.reshape(400, 1)))


def get_features(sequence):
    """
            Changing the rows and columns of amino acid frequency arrays

            Parameters
            -----------
            seq: str
            Sequence of phosphatase

            Returns
            ----------
            aac_array
            an array of amino acid frequencies in the form of 1 * 20
            """
    aac = getAAC(sequence)
    features = np.array(aac).reshape(1, 20)
    return features


def retrieve_model():
    """
    Train and Selection algorithm, and finally return an object containing the model object and
    the standardization tool (the sklearn standard scaler) object to return a tuple (ttopt, sscaler).
    """

    # Get dataset
    # filedata = os.path.join(this_dir, '../data/Topt249.tsv')
    # filedata = os.path.join(this_dir, '../data/Topt218.tsv')
    filedata = os.path.join(this_dir, '../data/Topt_rcaa218.tsv')
    with open(filedata, 'r', encoding='utf-8') as contrast:
        text = contrast.read()
        line = text.split('\n')
        contrast.close()
    for item in line:
        singledatalist.append(item.split('\t'))
    for i in singledatalist:
        alldatalist.append(i)

    # Get Dataset Information
    bar = tqdm(total=len(alldatalist), desc='Process completion progress1')
    for j in alldatalist[1:]:
        bar.update(1)
        namelist.append(j[2])
        seqlist.append(j[4])
        toptlist.append(j[6])

    aac = np.array([getAAC(seq) for seq in seqlist if seq != 'sequence'])
    aacdimer = np.array([get_dimer_frequency(seq) for seq in seqlist if seq != 'sequence'])

    for seq in seqlist:
        if seq != 'sequence':
            X = ProteinAnalysis(seq)
            instability_indexs.append(X.instability_index())
            molecular_weights.append(X.molecular_weight())
            aromaticities.append(X.aromaticity())
            secondaries.append(list(X.secondary_structure_fraction())[0])
            structures.append(list(X.secondary_structure_fraction())[1])
            fractions.append(list(X.secondary_structure_fraction())[2])

    # Change Feature Dimensions
    aromaticitiy = np.array(aromaticities)
    aromaticity = aromaticitiy[:, np.newaxis]
    instability_index = np.array(instability_indexs)
    instability_inde = instability_index[:, np.newaxis]
    molecular_weight = np.array(molecular_weights)
    molecularweight = molecular_weight[:, np.newaxis]
    secondariy = np.array(secondaries)
    secondariy = secondariy[:, np.newaxis]
    structure = np.array(structures)
    structure = structure[:, np.newaxis]
    fraction = np.array(fractions)
    fraction = fraction[:, np.newaxis]
    # Triad = np.array([list(GetProDes(seq).GetTriad().values()) for seq in seqlist if seq != 'sequence'])
    # PAAC = np.array([list(PseudoAAC.GetPseudoAAC(seq, lamda=100, weight=0.05, AAP=[_Hydrophobicity, _hydrophilicity]).values()) for seq in seqlist if seq != 'sequence'])  #
    # CTD = np.array([list(CalculateD(seq).values()) for seq in seqlist if seq != 'sequence'])
    stack1 = np.hstack((aac, molecularweight))
    bar.close()
    for topt in toptlist:
        if topt != 'temperature':
            labelslist.append(float(topt))

    sscaler = StandardScaler()

    # Processing Features and Labels
    X = sscaler.fit_transform(stack1)
    y = pd.Series(labelslist)

    # Select Algorithm and Strategy
    algorithm = KNeighborsRegressor(n_neighbors=1)
    Tpho = resreg.Rebagg(m=1, s=132, base_reg=algorithm)

    result = 0
    i = 1
    y_actual = []
    y_predicted = []
    for rrr in range(10):
        # Divide the test set and training set
        train_indices, test_indices = resreg.uniform_test_split(X, y, bins=bins, bin_test_size=0.25, random_state=rrr)
        X_train, y_train = X[train_indices, :], y[train_indices]
        X_test, y_test = X[test_indices, :], y[test_indices]
        # print(X_train, y_train)
        # print(X_test, y_test)
        final_feature = []
        final_feature.clear()
        '''
        ###
        for X_train_len in range(len(X_train)):
            feature_list = []
            feature_list.clear()
            for select_index in selected_index:
                # print('len_feature_list2')  # 400
                # print(len(feature_list))
                feature_list.append(X_train[X_train_len][select_index])
            # print('len_feature_list3')  # 400
            # print(len(feature_list))
            final_feature.append(feature_list)
        ###
        '''
        # print('final_feature_shape:', np.array((final_feature)).shape
        correlation = resreg.sigmoid_relevance(y_train, cl=21, ch=64)
        Tpho.fit(X_train, y_train, relevance=correlation, relevance_threshold=0.5,
                   sample_method='random_oversample', size_method='variation',
                   random_state=0)  # 重采样数据拟合回归器 （随机过采样）

        pred_test = Tpho.predict(X_test)
        # print('pred_test')
        # print(pred_test)
        y_actual += list(y_test)
        y_predicted += list(pred_test)
        result += sm.r2_score(y_test, pred_test)
        print(sm.r2_score(y_test, pred_test))
        i += 1
    print(f'result{result / 10}')

    # Draw a curve chart
    plt.plot(y_actual, label='Actual')
    plt.plot(y_predicted, label='Predicted')

    # Set coordinate axis range
    plt.xlim(1, 10)
    plt.ylim(0, 120)

    # Add Legend and Title
    plt.legend()
    plt.title('Actual vs Predicted')

    # Display graphics
    plt.show()

    return Tpho, sscaler


Tpho, sscaler = retrieve_model()


# Predict Topt of new sequences
# ================================#

def get_seq_name_fasta(fasta):
    """
    Obtain protein names and sequences in a fasta file
    
    Parameters
    -----------
    fasta: str
    Filename of fasta file containing protein names and sequences
    
    Returns
    ----------
    (names_list, sequences_list)
    A tuple of corresponding lists of  protein names
    and protein sequences in fasta_file
    """
    with open(fasta, 'r') as file:
        names, sequences = [], []
        for line in file:
            if line.startswith('>'):
                head = line.replace('>', '').strip()
                names.append(head)
                sequences.append('')
            else:
                seq = line.strip()
                if len(seq) > 0:
                    sequences[-1] += seq

    return names, sequences


def pred_phosphatase_seq_topt(sequence):
    """
    Predict the optimal catalytic temperature of phosphatase sequence.

    Parameters
    -----------
    sequence : str
        amino acid sequences of phosphatases

    Returns
    ---------
        Predicted optimal temperature for phosphatase catalysis.

    Examples
    ----------
    # >>> sequence = 'sequence'
    # >>> result = topt
    """
    moles = []

    # Calculate amino acid frequency and dipeptide frequency
    aac = np.array([getAAC(sequence)])

    # Obtaining protein information:
    X = ProteinAnalysis(sequence)
    moles.append(X.molecular_weight())

    # Change Feature Dimensions
    molee = np.array(moles)
    mole = molee[:, np.newaxis]

    # Obtain the final feature
    stack = np.hstack((aac, mole))

    # Standardized features
    X = sscaler.fit_transform(stack)

    y_pred = Tpho.predict(X)
    return y_pred


def pred_phosphatase_topt(fasta_file):

    """
    Predict the optimal catalytic temperature of phosphatase sequences in a fasta file.
    
    Parameters
    -----------
    fasta_file : str
        Fasta file containing amino acid sequences of phosphatases
    
    Returns
    ---------
        Predicted optimal temperature for phosphatase catalysis.
    
    Examples
    ----------
    # >>> fasta_file = 'test/sequences.fasta'
    # >>> result = pred_phosphatase_topt(fasta_file)
    """
    indexs = []
    moles = []
    structuress = []
    secondariess = []
    aromaticitiess = []
    fractionss = []
    names, sequences = get_seq_name_fasta(fasta_file)

    # Calculate amino acid frequency and dipeptide frequency
    aac = np.array([getAAC(sequence) for sequence in sequences])
    aacdimer = np.array([get_dimer_frequency(sequence) for sequence in sequences])

    # Obtaining protein information
    for sequence in sequences:
        X = ProteinAnalysis(sequence)
        indexs.append(X.instability_index())
        moles.append(X.molecular_weight())
        aromaticitiess.append(X.aromaticity())
        secondariess.append(list(X.secondary_structure_fraction())[0])
        structuress.append(list(X.secondary_structure_fraction())[1])
        fractionss.append(list(X.secondary_structure_fraction())[2])

    # Change Feature Dimensions
    molee = np.array(moles)
    mole = molee[:, np.newaxis]
    aromaticitiy = np.array(aromaticitiess)
    aromaticity = aromaticitiy[:, np.newaxis]
    secondariy = np.array(secondariess)
    secondariy = secondariy[:, np.newaxis]
    structure = np.array(structuress)
    structure = structure[:, np.newaxis]
    fraction = np.array(fractionss)
    fraction = fraction[:, np.newaxis]

    # Obtain the final feature
    stack = np.hstack((aac, mole))

    # Standardized features
    X = sscaler.fit_transform(stack)

    y_pred = Tpho.predict(X)
    pred_data = pd.DataFrame([names, list(y_pred)]).transpose()
    pred_data.columns = ['uniprot id', 'topt_pred']
    return pred_data


# if __name__ == "__main__":
#     # print(pred_phosphatase_seq_topt(sequence))
