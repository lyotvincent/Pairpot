import React from 'react';
import Markdown from 'react-markdown';
import ReactMarkdown from 'react-markdown';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { solarizedlight } from 'react-syntax-highlighter/dist/esm/styles/prism';

var Rcode = 
`
### Rcode
**********
\`\`\`r 
library(rhdf5)
library(SingleCellExperiment)
Load_h5adsc_to_SCE <- function(scmat, scgnm = NA){
  scmat <- h5read(scmat, '/')
  if("raw" %in% attr(scmat, "names")){
    X <- scmat$raw$X
    clsX <- class(X)
    if(length(clsX)>1)
      clsX <- clsX[1]
    if(clsX!="list"){
      X <- Matrix(X)
      dat <- as(X, "CsparseMatrix")
    }
    h5ad.var <- scmat$raw$var
  } else{
    X <- scmat$X
    dat <- sparseMatrix(i = X$indices[] + 1,
                        p = X$indptr[],
                        x = as.numeric(X$data[]),
                        repr = "C")
    h5ad.var <- scmat$var
  }
  # generate factors using categories
  var <- list()
  if("__categories" %in% attr(h5ad.var, "names")){ # old anndata
    for(name in attr(h5ad.var[["__categories"]], "names")){
      if(length(h5ad.var[[name]]) >= length(h5ad.var[["__categories"]][[name]])){
        var[[name]] <- factor(h5ad.var[[name]], labels = h5ad.var[["__categories"]][[name]])
      }
    }
  } else {
    for(name in attr(h5ad.var, "names")){
      if(name!='_index')
        if(class(var[[name]]) == "list"){
          var[[name]] <- factor(h5ad.var[[name]]$codes, labels = h5ad.var[[name]]$categories)
        }else{
          var[[name]] <- h5ad.var[[name]]
        }
    }
  }
  h5ad.obs <- scmat$obs
  obs <- list()
  if("__categories" %in% attr(h5ad.obs, "names")){ # old anndata
    for(name in attr(h5ad.obs[["__categories"]], "names")){
      if(length(h5ad.obs[[name]]) >= length(h5ad.obs[["__categories"]][[name]]))
        obs[[name]] <- factor(h5ad.obs[[name]], labels = h5ad.obs[["__categories"]][[name]])
    }
  } else{ # new anndata
    for(name in attr(h5ad.obs, "names")){
      if(name!='_index')
        if(class(h5ad.obs[[name]]) == "list"){
          obs[[name]] <- factor(h5ad.obs[[name]]$codes, labels = h5ad.obs[[name]]$categories)
        }else{
          obs[[name]] <- h5ad.obs[[name]]
        }
    }
  }

  dims <- c(length(h5ad.var[["_index"]]), length(h5ad.obs[["barcodes"]]))
  dat@Dim <- dims
  dat@Dimnames <- list(as.character(h5ad.var[["_index"]]), as.character(h5ad.obs[["_index"]]))
  obs <- data.frame(obs)
  rownames(obs) <- h5ad.obs[["barcodes"]] 
  return(sce)
}
\`\`\`
`

const sce = () => {
    return (
        <div style={{height:'25rem',overflow: 'auto'}}>
            {/* <ReactMarkdown>{Rcode}</ReactMarkdown> */}
            <ReactMarkdown
                children={Rcode}
                components={{
                code({ node, inline, className, children, ...props }) {
                    return inline ? (
                    <code className={className} {...props}>
                        {children}
                    </code>
                    ) : (
                    <SyntaxHighlighter
                        customStyle={{ backgroundColor: 'white' }}
                        style={solarizedlight}
                        language="r"
                        PreTag="div"
                        {...props}
                    >
                        {String(children)}
                    </SyntaxHighlighter>
                    );
                },
                }}
            />
        </div>
  );
};

export default sce;