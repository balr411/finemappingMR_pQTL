import scipy.sparse as sparse
import numpy as np
# Note that for indices, I can make this a file that it can read in 
def sparse_ld_subset(file_name, indices):
    data = data = np.load(file_name)
    idx = np.where(np.logical_and(np.isin(data['row'],indices),np.isin(data['col'],indices)))
    indices2seq = dict(zip(indices,range(len(indices))))
    irow = [indices2seq[x] for x in data['row'][idx]]
    icol = [indices2seq[x] for x in data['col'][idx]]
    R = sparse.coo_matrix((data['data'][idx], (irow, icol)), shape=(len(indices), len(indices))).toarray()
    R += R.T
    return R
