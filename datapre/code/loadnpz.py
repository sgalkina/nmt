from numpy import load

data = load('qm9_relgcn_kekulized_ggnp.npz')
lst = data.files
for item in lst:
    print(item)
    print(data[item])