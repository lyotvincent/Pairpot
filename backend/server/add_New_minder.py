from AddMENDER import AddMender4SPdata
import os
# path='/data/rzh/RawUrls'

path='/data/rzh/cncb/cncb'

for root, dirs, files in os.walk(path):
    for file in files:
        if file.startswith('New') and file.endswith('h5ad') and not 'sc' in root and not 'SC' in root and not 'sc' in file:
            print(os.path.join(root,file))
            filepath=os.path.join(root,file)
            AddMender4SPdata(filepath)