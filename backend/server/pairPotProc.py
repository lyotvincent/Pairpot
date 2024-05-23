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
        print("File ",directory[-3:]," is too large with {}".format(adata.shape[0]),file=file0)
        print("File ",directory[-3:]," is too large with {}".format(adata.shape[0]))
        return
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
    # adata = anno(adata, annoDict)
    # sc.pl.umap(adata, color='annotation')
    # adata = marker(adata)
    print(adata)
    adata.write_h5ad(directory+"/New_"+directory[-3:]+".h5ad")
    return

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
process_files('/data/rzh/RawUrls')
file0.close()
print("End")
f=open("End","w")
print("End",file=f)
f.close()