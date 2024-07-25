# evaluation the robustness of LPA
import h5py
import scipy.sparse
import time
import random
import numpy as np
import pandas as pd
import label_propagation as LPA
from sklearn.metrics.cluster import adjusted_rand_score

df = pd.DataFrame(columns=['MR', 'ARI_o', 'ARI_r', 'Time'])

start=time.time()
with h5py.File('resources/235/sp_deconv.h5ad','r') as f:
    group=f['obsp']['connectivities']

    data=group['data'][:]
    indices=group['indices'][:]
    indptr=group['indptr'][:]
    shape=(f['obsp']['connectivities'].attrs['shape'][0],f['obsp']['connectivities'].attrs['shape'][1])

    mat=scipy.sparse.csr_matrix((data,indices,indptr),shape=shape)


coo=mat.tocoo()

rows=coo.row
cols=coo.col
data=coo.data

with h5py.File('resources/235/sp_deconv.h5ad', 'r') as h5file:
    obs_group = h5file['obs']
    if "codes" in obs_group['annotation']:
      mat = obs_group['annotation']['codes'][:]
    else:
      mat = obs_group['annotation'][:]
val={}
for i in mat:
    if i not in val:
        val[i]=len(val)
        
k1=1.0
while True:
    for _ in range(100):
        X = LPA.matCoo(mat.shape[0], mat.shape[0])
        for i in range(len(data)):
            X.append(rows[i], cols[i], data[i])

        y_label = LPA.mat(mat.shape[0], len(val))
        y_new = LPA.mat(mat.shape[0], len(val))
        random_list=random.sample(range(mat.shape[0]), int(mat.shape[0] * 0.1))
        select_list=np.zeros(mat.shape[0])
        y_label.setneg()


        for i in random_list:
            select_list[i]=1
        for i in range(mat.shape[0]):
            if select_list[i]:
                y_label.editval2(i,val[mat[i]])
#case-ori
        start_time = time.perf_counter()
        
        y_pred = LPA.mat(mat.shape[0], len(val))
        y_res = LPA.mat(mat.shape[0], len(val))
        # LPA.dataProcess(y_label,y_new,k1,(1-k1)*0.5,(1-k1)*0.5)
        LPA.dataProcess(y_label,y_new,k1,(1-k1),0)
        LPA.labelPropagation(X, y_new, y_pred,y_res,0.5,1000)
        
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        
        res_arr = np.zeros(mat.shape[0])
        for i in range(mat.shape[0]):
            res_arr[i] = y_pred.getval(i,0)
        out_arr = np.array(mat)
        ari_o = adjusted_rand_score(out_arr, res_arr)
        # file2.write("{},{}\n".format((1-k1)*0.5,ari))
        # file2.write("{},{}\n".format((1-k1),ari_o))
        # print("{},{}".format((1-k1)*0.5,ari))
        
#case-rectified
        res_arr = np.zeros(mat.shape[0])
        for i in range(mat.shape[0]):
            res_arr[i] = y_res.getval(i,0)

        ari_r = adjusted_rand_score(out_arr, res_arr)
        # file3.write("{},{}\n".format((1-k1)*0.5,ari))
        # file3.write("{},{}\n".format((1-k1),ari_r))
        item = [round((1-k1),2), ari_o, ari_r, round(execution_time,5)]
        print(item)
        df.loc[len(df)] = item
        # print("{},{}".format((1-k1)*0.5,ari))
    k1-=0.05
    print(k1)
    if(k1<0):  break
df.to_csv("mistake-sp.tsv", sep='\t', header=True, index=True)
end=time.time()
print("time :{}".format(end-start))
