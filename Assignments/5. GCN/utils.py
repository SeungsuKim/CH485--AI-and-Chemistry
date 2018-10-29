import sys
from tqdm import tnrange, tqdm_notebook
import time

import numpy as np

from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP

from torch.utils.data import Dataset


def read_ZINC_smiles(file_name, num_mol):
    f = open(file_name, 'r')
    contents = f.readlines()

    smi_list = []
    logP_list = []

    for i in tqdm_notebook(range(num_mol), desc='Reading Data'):
        smi = contents[i].strip()
        m = Chem.MolFromSmiles(smi)
        smi_list.append(smi)
        logP_list.append(MolLogP(m))

    logP_list = np.asarray(logP_list).astype(float)

    return smi_list, logP_list


def smiles_to_onehot(smi_list):
    def smiles_to_vector(smiles, vocab, max_length):
        while len(smiles) < max_length:
            smiles += " "
        vector = [vocab.index(str(x)) for x in smiles]
        one_hot = np.zeros((len(vocab), max_length), dtype=int)
        for i, elm in enumerate(vector):
            one_hot[elm][i] = 1
        return one_hot

    vocab = np.load('./vocab.npy')
    smi_total = []

    for i, smi in tqdm_notebook(enumerate(smi_list), desc='Converting to One Hot'):
        smi_onehot = smiles_to_vector(smi, list(vocab), 120)
        smi_total.append(smi_onehot)

    return np.asarray(smi_total)


class OneHotLogPDataSet(Dataset):
    def __init__(self, list_one_hot, list_logP):
        self.list_one_hot = list_one_hot
        self.list_logP = list_logP

    def __len__(self):
        return len(self.list_one_hot)

    def __getitem__(self, index):
        return self.list_one_hot[index], self.list_logP[index]


def partition(list_one_hot, list_logP, args):
    num_total = list_one_hot.shape[0]
    num_train = int(num_total * (1 - args.test_size - args.val_size))
    num_val = int(num_total * args.val_size)
    num_test = int(num_total * args.test_size)

    one_hot_train = list_one_hot[:num_train]
    logP_train = list_logP[:num_train]
    one_hot_val = list_one_hot[num_train:num_train + num_val]
    logP_val = list_logP[num_train:num_train + num_val]
    one_hot_test = list_one_hot[num_total - num_test:]
    logP_test = list_logP[num_total - num_test:]

    train_set = OneHotLogPDataSet(one_hot_train, logP_train)
    val_set = OneHotLogPDataSet(one_hot_val, logP_val)
    test_set = OneHotLogPDataSet(one_hot_test, logP_test)

    partition = {
        'train': train_set,
        'val': val_set,
        'test': test_set
    }

    return partition