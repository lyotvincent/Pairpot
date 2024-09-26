import os
import zipfile

def tozip(file_path,zip_path):
    with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as zipf:
        zipf.write(file_path,os.path.basename(file_path))
        
rootpath="/home/rzh/stpair/backend/resources"
for root,dir,filelist in os.walk(rootpath):
    for filename in filelist:
        if filename.endswith('meta.h5ad'):
            file=os.path.join(root,filename)
            zip_file=file+".zip"
            print(f'file:{file},zip_file:{zip_file}')
            tozip(file,zip_file)