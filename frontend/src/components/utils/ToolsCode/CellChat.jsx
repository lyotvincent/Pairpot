import React from 'react';
import Markdown from 'react-markdown';
import ReactMarkdown from 'react-markdown';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { solarizedlight } from 'react-syntax-highlighter/dist/esm/styles/prism';

var Rcode = 
`
### CellChat
**Source:** https://github.com/sqjin/CellChat\n
**Citation:** Jin, S., Guerrero-Juarez, C.F., Zhang, L., Chang, I., Ramos, R., Kuan, C., Myung, P., Plikus, M.V., and Nie, Q. (2021). Inference and analysis of cell-cell communication using CellChat. Nat. Commun. 12, 1088. 
### Key functions of CellChat include:
- Identifying Ligand-Receptor Interactions
- Inferring Cell Communication Networks 

### Rcode
**********
\`\`\`r 
# fetching params
library(getopt)
spec <- matrix(c(
  "help", "h", 0, "logical", "help",
  "path", "p", 1, "character", "path to spatial or single-cell h5ad file",
  "outpath", "o", 1, "character", "path to output file, this script will generate two dataframes",
  "species", "s", 2, "character", "species for CellChat, Mouse or Human, default Mouse.",
  "raw", "r", 2, "logical", "whether use raw data in h5ad file"
), byrow = T, ncol = 5)
args <- getopt(spec)
if (!is.null(args$help) || is.null(args$path) || is.null(args$outpath)) {
  cat(paste(getopt(spec, usage = T), "\\n"))
  q(status = 1)
}

species <- args$species
if (is.null(species)) {
  species <- "Mouse"
}

raw <- args$raw
if (is.null(raw)) {
  raw <- TRUE
}

source("_loading.R")
library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
library(stringr)
library(ggplot2)
library(NMF)

LegRec_CellChat <- function(spe, species = "Mouse", is_cell = FALSE) {
  sp_data <- spe@assays@data@listData$counts
  colnames(sp_data) <- Make_unique(colnames(sp_data))
  rownames(sp_data) <- str_to_title(rownames(sp_data))
  sp_celltype <- data.frame(
    annotation = spe$annotation,
    row.names = colnames(sp_data)
  )

  # filter the 'NA' labels
  sp_data <- sp_data[, which(!is.na(sp_celltype$annotation))]
  sp_celltype <- data.frame(
    annotation = sp_celltype[which(!is.na(sp_celltype$annotation)), ],
    row.names = colnames(sp_data)
  )

  # set the labels from 0+ to 1+
  do_add <- FALSE
  if ("0" %in% levels(sp_celltype$annotation)) {
    levels(sp_celltype$annotation) <- as.character(as.numeric(levels(sp_celltype$annotation)) + 1)
    do_add <- TRUE
  }

  # create cellchat object
  cellchat <- createCellChat(object = sp_data,
                             meta = sp_celltype, 
                             group.by = "annotation")
  groupSize <- as.numeric(table(cellchat@idents))
  species <- "Mouse"
  # Load the cellchatDB
  if (species == "Mouse") {
    CellChatDB <- CellChatDB.mouse
  } else if (species == "Human") {
    CellChatDB <- CellChatDB.human
  }

  # The corresponding classification is taken out and used as the analysis database
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  # Load the database to 'cellchat' object
  cellchat@DB <- CellChatDB
  # Extract expression data
  cellchat <- subsetData(cellchat)
  # Search hyper variable genes
  cellchat <- identifyOverExpressedGenes(cellchat)
  # Search hyper pathways
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # mapping data to PPI and save the results to 'cellchat@LR$LRsig'
  cellchat <- projectData(cellchat, PPI.human)
  # Calculate the probability of cell-cell communication
  cellchat <- computeCommunProb(cellchat, raw.use = T)

  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  df <- subsetCommunication(cellchat)
  net <- cellchat@net$count
  if (do_add) {
    df["source"] <- df["source"] - 1
    df["target"] <- df["target"] - 1
    rownames(net) <- as.character(as.numeric(rownames(net)) - 1)
    colnames(net) <- as.character(as.numeric(colnames(net)) - 1)
  }
  output <- list(heatmap = df, network = net)
  return(output)
}

Run_CellChat <- function(inpath, outpath, species, is_cell = F, raw = T) {
  cat("**Loading h5ad data to SingleCellExperiment...\\n")
  sce <- Load_h5adsc_to_SCE(inpath, raw = raw)
  print(sce)
  cat("**Loading data finished, running CellChat...\\n")
  output <- LegRec_CellChat(sce, species, is_cell)
  cat("**Cellchat finished, writing to dictionary...\\n")
  write.table(output$heatmap,
              paste0(outpath, "/CellChat_heatmap.tsv"),
              sep = "\\t",
              quote = FALSE)
  write.table(output$net,
              paste0(outpath, "/CellChat_network.tsv"),
              sep = "\\t",
              quote = FALSE)
  message(paste0(
    "CCI with 'CellChat' finished, idents saved at ",
    paste0(outpath, "/CellChat_heatmap.tsv"), " and ",
    paste0(outpath, "/CellChat_network.tsv")
  ))
}

cat("#### Running CellChat ####\\n")
cat(paste0("**Input h5ad file:", args$path, "\\n"))
cat(paste0("**Output h5ad file:", args$outpath, "\\n"))
cat(paste0("**Species:", species, "\\n"))
cat(paste0("**Using raw:", as.character(raw), "\\n"))
Run_CellChat(args$path, args$outpath, species, is_cell = FALSE, raw = raw)
cat("#### CellChat Finished ####\\n")
\`\`\`
`

const CellChat = () => {
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

export default CellChat;