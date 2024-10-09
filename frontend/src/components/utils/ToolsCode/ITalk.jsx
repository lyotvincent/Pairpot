import React from 'react';
import Markdown from 'react-markdown';
import ReactMarkdown from 'react-markdown';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { solarizedlight } from 'react-syntax-highlighter/dist/esm/styles/prism';

var Rcode = 
`
### iTALK
iTALK is a computational tool designed to analyze and visualize intercellular communication, particularly focusing on interactions between cells through ligand-receptor pairs.

### Key features and functions of iTALK include:
- Identifying Ligand-Receptor Interactions
- Visualizing Cell-Cell Communication

### Rcode
**********
\`\`\`r 
# fetching params
library(getopt)
spec <- matrix(c(
  "help", "h", 0, "logical", "help",
  "path", "p", 1, "character", "path to spatial or single-cell h5ad file",
  "outpath", "o", 1, "character", "path to output file, 
  this script will generate two dataframes",
  "raw", "r", 2, "logical", "whether use raw data in h5ad file"
), byrow = TRUE, ncol = 5)
args <- getopt(spec)
if (!is.null(args$help) || is.null(args$path) || is.null(args$outpath)) {
  cat(paste(getopt(spec, usage = TRUE), "\\n"))
  q(status = 1)
}

raw <- args$raw
if (is.null(raw)) {
  raw <- TRUE
}

source("_loading.R")
library(iTALK)

LegRec_iTALK <- function(sce) {
  exp <- as.data.frame(t(as.matrix(sce@assays@data$counts)))
  exp$cell_type <- as.character(sce$annotation)
  highly_exprs_genes <- rawParse(exp, top_genes = 50, stats = "mean")
  comm_list <- c("growth factor", "other", "cytokine", "checkpoint")
  res <- NULL
  for (comm_type in comm_list) {
    res_cat <- FindLR(highly_exprs_genes,
                      datatype = "mean count",
                      comm_type = comm_type)
    res_cat <- res_cat[which(res_cat$cell_from_mean_exprs *
                               res_cat$cell_to_mean_exprs > 0), ]
    res_cat <- res_cat[order(res_cat$cell_from_mean_exprs *
                               res_cat$cell_to_mean_exprs,
                             decreasing = TRUE), ]
    res <- rbind(res, res_cat)
  }
  output <- list(heatmap = res)
  return(output)
}

Run_iTALK <- function(inpath, outpath, raw) {
  sce <- Load_h5adsc_to_SCE(inpath, raw = raw)
  rownames(sce) <- toupper(rownames(sce))
  rownames(sce@assays@data@listData$counts) <- toupper(rownames(sce@assays@data@listData$counts))
  output <- LegRec_iTALK(sce)
  write.table(output$heatmap,
              paste0(outpath, "/iTALK_heatmap.tsv"),
              sep = "\\t",
              quote = FALSE)
}

cat("#### Running iTALK ####\\n")
cat(paste0("**Input h5ad file:", args$path, "\\n"))
cat(paste0("**Output h5ad file:", args$outpath, "\\n"))
cat(paste0("**Using raw:", as.character(raw), "\\n"))
Run_iTALK(args$path, args$outpath, raw = raw)
cat("#### iTALK Finished ####\\n")
\`\`\`
`

const iTALK = () => {
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

export default iTALK;