# https://bioconductor.org/books/3.18/OSCA.workflows/zeisel-mouse-brain-strt-seq.html

stop("You don't need to run this for the workshop.") 
    
library(scRNAseq)
library(scater)
library(org.Mm.eg.db)
library(scran)

sce = ZeiselBrainData()

sce <- aggregateAcrossFeatures(sce, 
                               id=sub("_loc[0-9]+$", "", rownames(sce)))

rowData(sce)$Ensembl <- mapIds(org.Mm.eg.db, 
                               keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")

unfiltered <- sce

stats <- perCellQCMetrics(sce, subsets=list(
    Mt=rowData(sce)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", 
                                              "subsets_Mt_percent"))
sce <- sce[,!qc$discard]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

set.seed(1000)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters) 
sce <- logNormCounts(sce)

dec.zeisel <- modelGeneVarWithSpikes(sce, "ERCC")

top.hvgs <- getTopHVGs(dec.zeisel, n = 500)

# add fake dropouts in one cell type --------------------------------------

library(tinyplot); tinytheme("clean2", family = "Arial")
library(fastverse)

options(digits = 3)

m = sce[top.hvgs,] |> 
    logcounts() |> 
    t()

ur = uwot::umap(m) |> 
    set_colnames(c("umap_1", "umap_2")) |> 
    qDT(row.names.col = "cell_id") |> 
    mtt(l1c = sce$level1class |> forcats::fct())

plt(umap_2 ~ umap_1 | l1c,
    data = ur, 
    pch = 15,
    cex=.2,
    legend = legend(pt.cex=1),
    main = "Original data UMAP")

dat = ur |> 
    add_vars(m |> qDT()) 

(g2d = dat |> # genes to dropout
    pivot(ids = 1:4) |> 
    sbt(l1c %like% "CA1") |> 
    sbt(!(variable %like% "OTTMUSG|mt-")) |> 
    gby(variable) |> 
    smr(var = fvar(value)) |> 
    roworder(var) |> 
        tail(20))

(to_imp = as.character(g2d$variable))

# In each of the 10 most highly expressed genes in pyramidal CA1 cells, set 20% of observations to NA
add_NAs = function(x, ct) {
    x[ct %==% "pyramidal CA1"][sample.int(938, 188)] = NA
    x
}

set.seed(123)

na_dat = dat |> 
    mtt(across(to_imp,
               add_NAs,
               ct = l1c))

na_dat |> na.omit() |> dim()

na_dat |> na.omit() |> fcount(l1c) # few CA1 cells left

na_umap = na_dat |> 
    slt(-(cell_id:l1c)) |> 
    qM() |> 
    na.omit() |> 
    uwot::umap() |> 
    qDT() |> 
    set_colnames(c("umap_1", "umap_2")) |> 
    mtt(na_dat |> na.omit() |> slt(l1c)) 

na_dat = na_dat |> 
    slt(-(umap_1:umap_2))
    
plt(umap_2 ~ umap_1 | l1c,
    data = na_umap,
    pch = 15,
    cex = .33,
    legend = legend(pt.cex=1),
    main = "NAs omitted")

imp_input = na_dat |> 
    qDF() |> 
    slt(-(cell_id:l1c)) |> 
    setRownames(na_dat$cell_id) 

head(imp_input[,1:5])

imp = mice(m = 1,
           imp_input, 
           method = "norm")

imp_dat = copy(imp_input)

# Copy the imputated values into their respective places in imp_dat:
for (g in to_imp) {
    g_imp = imp$imp[[g]]
    
    imp_dat[[g]][fmatch(rownames(g_imp), rownames(imp_dat))] = g_imp$`1`
}

imp_dat |> slt(1:5)

imp_dat |> na.omit() |> nrow() # none dropped now

set.seed(123)

imp_umap = imp_dat |> 
    qM() |> 
    uwot::umap() |> 
    qDT() |> 
    set_colnames(c("umap_1", "umap_2")) |> 
    mtt(na_dat |> slt(l1c)) 

plt(umap_2 ~ umap_1 | l1c,
    data = imp_umap,
    pch = 15,
    cex = .33,
    legend = legend(pt.cex=1),
    main = "Imputation seems to recover missing values\nin pyramidal CA1 cells")

save(dat, na_dat, to_imp,
     file = "~/projects/missing_data/workshop/sc.RData",
     compress = "xz",
     compression_level = 9)
