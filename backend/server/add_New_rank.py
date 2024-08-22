from AddRank import AddRank4SPdata
import os
# path='/data/rzh/RawUrls'
path='/data/rzh/RawUrls/209/STDS0000209'
errors=[]
for root,dirs,files in os.walk(path):
    if not 'sc' in root and not 'SC' in root and not 'sc' in files:
        clu_key = "annotation"
    else:
        clu_key = "leiden-1"
    if "meta.txt" in files:
        print("root:",root)
        lines=[]
        metapath=os.path.join(root,"meta.txt")
        with open(metapath) as f:
            lines=f.readlines()
        lines=[line.strip() for line in lines]
        print(lines)
        for file in files:
            if file.startswith('New') and file.endswith('h5ad'):
                print(file)
                filepath=os.path.join(root,file)
                # try:
                #     AddRank4SPdata(spH5adFile=filepath,organs=lines, clu_key=clu_key)
                # except:
                #     errors.append(os.path.join(root,file))
                AddRank4SPdata(spH5adFile=filepath,organs=lines, clu_key=clu_key)
                
            # print("file:",file)
print("\n\n\n Error:")
for _ in errors:
    print(_)