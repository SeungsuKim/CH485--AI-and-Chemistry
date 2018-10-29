import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from torch.utils.data import Dataset
from sklearn.model_selection import train_test_split, StratifiedKFold

class ToxDataset(Dataset):
    def __init__(self, fps, toxs):
        self.fps = fps
        self.toxs = toxs
        
    def __len__(self):
        return len(self.fps)
    
    def __getitem__(self, index):
        return self.fps[index], self.toxs[index]
    

def read_data(filename):
    f = open(filename + '.smiles', 'r')
    contents = f.readlines()

    smiles = []
    labels = []
    for i in contents:
        smi = i.split()[0]
        label = int(i.split()[2].strip())

        smiles.append(smi)
        labels.append(label)

    num_total = len(smiles)
    rand_int = np.random.randint(num_total, size=(num_total,))
    
    return np.asarray(smiles)[rand_int], np.asarray(labels)[rand_int]

def get_fingerprint(smile, args):
    mol = Chem.MolFromSmiles(smile)
    arr = np.zeros((1,))
    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, args.n_bits)
        DataStructs.ConvertToNumpyArray(fp, arr)
    except:
        pass
    
    return arr

def make_dataset(args):
    smiles, labels = read_data(args.tox_type)
    
    fps = list()
    toxs = list()
    for i, smile in enumerate(smiles):
        #try:
            fp = get_fingerprint(smile, args)
            fps.append(fp)
            toxs.append(labels[i])
        #except:
        #    pass
    
    return np.array(fps), np.array(toxs)

def partition(fps, toxs, args):
    num_total = fps.shape[0]
    num_train = int(num_total*(1-args.test_size-args.val_size))
    num_val = int(num_total*args.val_size)
    num_test = int(num_total*args.test_size)
    
    fps_train = fps[:num_train]
    toxs_train = toxs[:num_train]
    fps_val = fps[num_train:num_train+num_val]
    toxs_val = toxs[num_train:num_train+num_val]
    fps_test = fps[num_total-num_test:]
    toxs_test = toxs[num_total-num_test:]
        
    train_set = ToxDataset(fps_train, toxs_train)
    val_set = ToxDataset(fps_val, toxs_val)
    test_set = ToxDataset(fps_test, toxs_test)
        
    partition = {
        'train': train_set,
        'val': val_set,
        'test': test_set
    }
    
    return partition
'''
def partition(fps, toxs, args):
    folded_partition = list()
    sfk = StratifiedKFold(n_splits=args.n_splits, random_state=arg.seed)
    for trains, tests in sfk.split(fps, toxs):
'''        