# Each process fits one method and one PC count on frozen matrices.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
pc <- as.integer(args[1])
prepared <- args[2]
out <- args[3]
stopifnot(pc %in% c(3L,5L,10L))
lib <- Sys.getenv("SV_R_LIBRARY")
if (nzchar(lib)) .libPaths(c(lib,.libPaths()))
cfg <- jsonlite::fromJSON(file.path(prepared,"settings.json"))
write_tsv <- function(data, path) {
    p <- file.path(out,path)
    connection <- if (endsWith(path,".gz")) gzfile(p,"wt") else file(p,"wt")
    on.exit(close(connection))
    write.table(data,connection,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
read_frozen <- function(prefix, matrix_name) {
    classes <- if (cfg$rda_blocked) c(id="character",block="character") else c(id="character")
    meta <- read.delim(file.path(prepared,paste0(prefix,"_design.tsv")),row.names=1,check.names=FALSE,
                      colClasses=classes)
    y <- as.matrix(read.delim(file.path(prepared,matrix_name),row.names=1,check.names=FALSE,
                             colClasses=c(id="character")))
    stopifnot(identical(rownames(y),rownames(meta)),!anyDuplicated(colnames(y)),all(is.finite(y)))
    list(y=y,meta=meta)
}
env_names <- cfg$environments
pc_names <- paste0("PC",seq_len(pc))
geo <- if (cfg$adjust_geography) c("lon","lat") else character(0)
cov_names <- c(pc_names,geo)
design <- function(meta) {
    x <- scale(as.matrix(meta[,c(env_names,cov_names),drop=FALSE]))
    if (!all(is.finite(x)) || qr(cbind(1,x))$rank!=1+ncol(x))
        stop("Design contains a constant or redundant environment/covariate.")
    x
}

d <- read_frozen("group","rda_frequency.tsv.gz")
x <- as.data.frame(design(d$meta))
Y <- d$y
form <- as.formula(paste("Y ~",paste(env_names,collapse=" + "),
                         "+ Condition(",paste(cov_names,collapse=" + "),")"))
fit <- vegan::rda(form,data=x,scale=TRUE)
df <- nrow(x)-1-ncol(x)
stopifnot(df>0,fit$CCA$rank>0L)
set.seed(cfg$seed)
restriction <- if (cfg$rda_blocked) factor(d$meta$block) else NULL
control <- permute::how(nperm=cfg$permutations,blocks=restriction)
permutation_matrix <- permute::shuffleSet(nrow(Y),nset=cfg$permutations,control=control)
if (nrow(permutation_matrix)<19L) stop("Too few valid group permutations. Review the block design.")
if (cfg$rda_blocked) {
    for (i in seq_len(nrow(permutation_matrix)))
        stopifnot(all(restriction[permutation_matrix[i,]]==restriction))
}
global <- anova(fit,model="reduced",permutations=permutation_matrix,parallel=1)
marginal <- anova(fit,by="margin",model="reduced",permutations=permutation_matrix,parallel=1)
by_axis <- anova(fit,by="axis",model="reduced",permutations=permutation_matrix,parallel=1)
axes <- names(fit$CCA$eig)
p_axis <- by_axis[match(axes,rownames(by_axis)),"Pr(>F)"]
q_axis <- p.adjust(p_axis,"BH")
p_global <- global[1,"Pr(>F)"]
eligible <- is.finite(q_axis)&q_axis<cfg$alpha&p_global<cfg$alpha
r2 <- vegan::RsquareAdj(fit)
write_tsv(data.frame(pc_count=pc,groups=nrow(Y),loci=ncol(Y),residual_df=df,
    permutations=nrow(permutation_matrix),blocked=cfg$rda_blocked,global_p=p_global,
    pseudo_f=global[1,"F"],constrained_fraction=fit$CCA$tot.chi/fit$tot.chi,
    adjusted_r2=r2$adj.r.squared),"global.tsv")
write_tsv(data.frame(axis=axes,inertia=as.numeric(fit$CCA$eig),p=p_axis,q=q_axis,eligible=eligible),
          "axes.tsv")
write_tsv(data.frame(environment=env_names,p=marginal[match(env_names,rownames(marginal)),"Pr(>F)"]),
          "marginal_terms.tsv")
choices <- seq_along(axes)
scores <- function(which) as.matrix(vegan::scores(fit,display=which,choices=choices,scaling=2))
loadings <- scores("species")
loci <- rownames(loadings)
z <- scale(loadings)
peak <- if (any(eligible)) apply(abs(z[,eligible,drop=FALSE]),1,max) else rep(NA_real_,length(loci))
result <- data.frame(locus=loci,max_abs_loading_z=peak,candidate=is.finite(peak)&peak>cfg$loading_z)
write_tsv(cbind(result,as.data.frame(loadings)),"locus_loadings.tsv.gz")
write_tsv(result[result$candidate,],"exploratory_candidates.tsv")
si <- scores("sites")
bp <- scores("bp")
write_tsv(data.frame(group=rownames(si),si),"site_scores.tsv")
write_tsv(data.frame(term=rownames(bp),bp),"environment_arrows.tsv")
write_tsv(data.frame(permutation_matrix),"permutations.tsv.gz")
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))
