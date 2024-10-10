from .lassoview import LPARefine
from .pairview import NNLSDeconvolution
from .pairPotProc import process_files
from .normalize import rank
from .download import downLd
import pandas as pd
import anndata as ad


def lassoView(selected,adata,use_model,do_correct=True,function="anndata"):
    result=LPARefine(selected=selected,adata=adata,use_model=use_model,do_correct=do_correct,function=function)
    return result

def lassoView(selected,adata,do_correct=True,function="anndata"):
    result=LPARefine(selected=selected,adata=adata,do_correct=do_correct,function=function)
    return result

def pairView(selected,scdata,spdata,function="anndata"):
    result=NNLSDeconvolution(selected=selected, scdata=scdata, spdata=spdata,function=function)
    return result

def lassoProc(path):
    result=process_files(path=path)
    return result

def pairProc(adata:ad.AnnData, organs, top=0.05, alpha=10e-40, n_jobs=16, clu_key='leiden-1'):
    result=rank(adata=adata,organs=organs,top=top,alpha=alpha,n_jobs=n_jobs,clu_key=clu_key)
    return result

def downLd(dataset_id,type,file):
    result=downLd(dataset_id=dataset_id,type=type,file=file)
    return dataset_id