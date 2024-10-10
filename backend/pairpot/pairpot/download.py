import wget
import anndata as ad
import zipfile
import os
import shutil

def extract_zip(zip_file_path, extract_to_path):
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to_path)
        
def downLd(dataset_id,type,file):
    bioxai=""
    if len(dataset_id>=3):
        dataset_id=dataset_id[-3:]
    bioxai="http://src.bioxai.cn/000"[:-len(dataset_id)]+dataset_id+"/"
    url=bioxai
    filename=""
    if type=="sc":
        if file=="meta":
            filename="sc_meta.h5ad.zip"
        elif file=="complete":
            filename="sc_sampled.h5ad"
        else:
            print('Wrong file, use "meta" or "complete" instead.')
            raise TypeError
    elif type=="sp":
        if file=="meta":
            filename="sp_meta.h5ad.zip"
        elif file=="complete":
            filename="sp_deconv.h5ad"
        else:
            print('Wrong file, use "meta" or "complete" instead.')
            raise TypeError
    else:
        print('Wrong type, use "sc" or "sp" instead.')
        raise TypeError
    print(url+filename)
    downld_url=url+filename
    
    path="./pairpot-download/"
    if "pairpot-download" not in os.listdir("./"):
        os.mkdir(path)
    else:
        shutil.rmtree(path)
        os.mkdir(path)
    wget.download(downld_url,path)
    if file=="meta":
        extract_zip(path+filename,path+filename[:-4])
        return ad.read_h5ad(path+filename[:-4]+"/"+filename[:-4]) 
    elif file=="complete":
        return ad.read_h5ad(path+filename)    
    else:
        print('Wrong file, use "meta" or "complete" instead.')
        raise TypeError