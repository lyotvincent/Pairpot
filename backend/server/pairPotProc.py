#nohup python '/home/rzh/stpair/backend/server/pairPotProc.py' >'/home/rzh/stpair/backend/server/ProcLog.txt' 2>&1 &
#source ~/Browser/bin/activate
#cd '/home/rzh/stpair/backend/server'

from normalize import *
import os
from glob import glob

file0=open("Error.txt","w")

def proc_h5ad(filename, filepath):
    directory = os.path.dirname(filepath[0])
    meta_file_path = os.path.join(directory, 'meta.txt')
    organs=[]
    if os.path.exists(meta_file_path):
        with open(meta_file_path, 'r') as file:
            organs = file.read().splitlines()
    print(filename,filepath,organs)
    print(directory+"/New_"+directory[-3:]+".h5ad")
    # samples=filepath
    # sampleNames=filename
    adata=concat_adata(filepath, filename, inputFunc=input_adata_h5ad)
    adata = pp(adata)
    if adata.shape[0]>100000:
        sampleIdx = np.random.choice(len(adata), 50000)
        adata = adata[sampleIdx,:].copy()
    adata = clu(adata)
    
    try:
        adata = rank(adata, organs)
    except:
        try:
            print("Try kruskal in",filepath[0])
            adata=rank(adata,organs,test_func=kruskal)
        except:
            print("Error in",filepath[0],file=file0)
            print("Error in",filepath[0])
    # print(adata.uns['UCell_Assign'])
    adata = marker(adata, groupby="leiden-1", method='wilcoxon')
    adata.obs = adata.obs.astype('str').fillna("")
    # adata = anno(adata, annoDict)
    # sc.pl.umap(adata, color='annotation')
    # adata = marker(adata)
    print(adata)
    adata.write_h5ad(directory+"/New_"+directory[-3:]+".h5ad")
    return

def proc_txtgz(filename, filepath, type='sc'):
    directory = os.path.dirname(filepath[0])
    meta_file_path = os.path.join(directory, 'meta.txt')
    organs=[]
    if os.path.exists(meta_file_path):
        with open(meta_file_path, 'r') as file:
            organs = file.read().splitlines()
    print(filename,filepath,organs)
    print(f"{directory}/{type}{directory[-3:]}_raw.h5ad")
    samples=filepath
    sampleNames=filename
    adata = concat_adata(samples, sampleNames, inputFunc=input_adata_txt)
    adata = pp(adata)
    adata = clu(adata)
    if adata.shape[0]>100000:
        print("File ",directory[-3:]," is too large with {}".format(adata.shape[0]),file=file0)
        print("File ",directory[-3:]," is too large with {}".format(adata.shape[0]))
        return
    try:
        adata = rank(adata, organs)
    except:
        try:
            print("Try kruskal in",filepath[0])
            adata=rank(adata,organs,test_func=kruskal)
        except:
            print("Error in",filepath[0],file=file0)
            print("Error in",filepath[0])
    adata = marker(adata, groupby="leiden-1", method='wilcoxon')
    adata.write_h5ad(f"{directory}/{type}{directory[-3:]}_raw.h5ad")

def proc_h5(filename, filepath):
    directory = os.path.dirname(filepath[0])
    meta_file_path = os.path.join(directory, 'meta.txt')
    organs=[]
    if os.path.exists(meta_file_path):
        with open(meta_file_path, 'r') as file:
            organs = file.read().splitlines()
    print(filename,filepath,organs)
    print(directory+"/New_"+directory[-3:]+".h5ad")
    samples=filepath
    sampleNames=filename
    adata = concat_adata(samples, sampleNames, inputFunc=input_adata_10Xh5)
    adata = pp(adata)
    adata = clu(adata)
    if adata.shape[0]>100000:
        print("File ",directory[-3:]," is too large with {}".format(adata.shape[0]),file=file0)
        print("File ",directory[-3:]," is too large with {}".format(adata.shape[0]))
        return
    try:
        adata = rank(adata, organs)
    except:
        try:
            print("Try kruskal in",filepath[0])
            adata=rank(adata,organs,test_func=kruskal)
        except:
            print("Error in",filepath[0],file=file0)
            print("Error in",filepath[0])
    # print(adata.uns['UCell_Assign'])
    adata = marker(adata, groupby="leiden-1", method='wilcoxon')
    # adata = anno(adata, annoDict)
    # sc.pl.umap(adata, color='annotation')
    # adata = marker(adata)
    # print(adata)
    adata.write_h5ad(directory+"/New_"+directory[-3:]+".h5ad")
    return


def process_files(path):
    for root, dirs, files in os.walk(path):
        # print(files)
        filepath=[]
        filename=[]
        flag=0
        for file in files:
            if(file.startswith('New_')):
                print(file)
                flag=1
                continue
            if(not(file.endswith('.h5ad') or file.endswith('h5'))):
                print(file)
                continue
            filepath.append(os.path.join(root,file))
            filename.append(file)
        # print(filepath)
        # print(filename)
        if filename==[] :
            continue
        if flag==1:
            continue
        if filename[0].endswith('.h5ad'):
            try:
                proc_h5ad(filename,filepath)
            except:
                print("Error in",filepath[0],file=file0)
                print("Error in",filepath[0])
        elif filename[0].endswith('.h5'):
            try:
                proc_h5(filename,filepath)
            except:
                print("Error in",filepath[0],file=file0)
                print("Error in",filepath[0])
        # for file in files:
        #     if file.endswith('.h5ad'):
        #         filepath = os.path.join(root, file)
                # proc_h5ad(file, filepath)
                
                # try:
                #     proc_h5ad(file, filepath)
                # except:
                #     print("Error in "+filepath)
            # elif file.endswith('.h5'):
            #     filepath = os.path.join(root, file)
                # proc_h5(file, filepath)
                # try:
                    # proc_h5(file, filepath)
                # except:
                #     print("Error in "+filepath)
print(1)
# process_files('/data/rzh/RawUrls')
# sample = {
#    'dictionary':'/data/rzh/RawUrls/227/SCDS0000005',
#    'path':["GSM4089151_P1_gene_cell_exprs_table.h5ad",
#           "GSM4089152_P2_gene_cell_exprs_table.h5ad",
#           "GSM4089153_P3_gene_cell_exprs_table.h5ad",
#           "GSM4089154_P4_gene_cell_exprs_table.h5ad",
#           "GSM4711414_P5_gene_cell_exprs_table.h5ad",
#           "GSM4711415_P6_gene_cell_exprs_table.h5ad"],
#     'name':['P1', 'P2', 'P3', 'P4', 'P5', 'P6']
# }
import json
sample = {
  "dict":"/data/rzh/RawUrls/212/SCDS0000002",
  "path":["GSM4653855_AD1.h5ad",
          "GSM4653856_AD2.h5ad",
          "GSM4653857_AD3.h5ad",
          "GSM4653858_AD4.h5ad",
          "GSM4653859_AD5.h5ad",
          "GSM4653860_AD6.h5ad",
          "GSM4653861_AD7.h5ad",
          'GSM4653862_AD8.h5ad',
          'GSM4653863_HC1.h5ad',
          'GSM4653864_HC2.h5ad',
          'GSM4653865_HC3.h5ad',
          'GSM4653866_HC4.h5ad',
          'GSM4653867_HC5.h5ad',
          'GSM4653868_HC6.h5ad',
          'GSM4653869_HC7.h5ad',],
  "name":["AD1",
          "AD2",
          "AD3",
          "AD4",
          "AD5",
          "AD6",
          "AD7",
          'AD8',
          'HC1',
          'HC2',
          'HC3',
          'HC4',
          'HC5',
          'HC6',
          'HC7',]
}
json_data = json.dumps(sample)
 
# 将JSON数据保存到文件
with open(f"{sample['dict']}/sample.json", 'w') as f:
    json.dump(sample, f)
proc_h5ad(filepath=[f"{sample['dict']}/{s}" for s in sample["path"]], filename=sample['name'])
file0.close()
print("End")
f=open("End","w")
print("End",file=f)
f.close()