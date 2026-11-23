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

d <- read_frozen("individual","lfmm_dosage.tsv.gz")
x <- design(d$meta)
g <- d$y
fit <- LEA::lfmm2(input=g,env=x,K=cfg$lfmm_k,lambda=1e-5,effect.sizes=FALSE)
model <- cbind(Intercept=1,x,fit@U)
qr_model <- qr(model)
df <- nrow(model)-ncol(model)
if (qr_model$rank != ncol(model) || df <= 0) stop("LFMM2 test design is rank deficient or saturated.")
stderr_unit <- numeric(ncol(model))
stderr_unit[qr_model$pivot] <- sqrt(diag(chol2inv(qr.R(qr_model))))
m <- ncol(g)
indices <- seq_along(env_names)+1L
coefficients <- statistics <- matrix(NA_real_,length(env_names),m)
for (first in seq.int(1L,m,by=2048L)) {
    j <- first:min(first+2047L,m)
    y <- g[,j,drop=FALSE]
    beta <- qr.coef(qr_model,y)
    rss <- colSums(qr.resid(qr_model,y)^2)
    tss <- colSums(sweep(y,2,colMeans(y))^2)
    ok <- is.finite(rss)&rss>1e-12*tss
    if (any(ok)) {
        coefficients[,j[ok]] <- beta[indices,ok,drop=FALSE]
        errors <- outer(stderr_unit[indices],sqrt(rss[ok]/df))
        statistics[,j[ok]] <- beta[indices,ok,drop=FALSE]/errors
    }
}
results <- diagnostics <- vector("list",length(env_names))
for (i in seq_along(env_names)) {
    stat <- statistics[i,]
    good <- is.finite(stat)
    inflation <- if (any(good)) median(stat[good]^2)/qchisq(.5,1) else NA_real_
    p_gc <- rep(NA_real_,m)
    reasons <- ifelse(good,"tested","degenerate_residual")
    if (is.finite(inflation) && inflation>0) {
        p_gc[good] <- pchisq(stat[good]^2/inflation,1,lower.tail=FALSE)
    } else {
        reasons[good] <- "invalid_genomic_control"
    }
    q <- p.adjust(p_gc,"BH")
    results[[i]] <- data.frame(locus=colnames(g),environment=env_names[i],pc_count=pc,
        latent_k=cfg$lfmm_k,beta=coefficients[i,],statistic=stat,p_ols=2*pt(-abs(stat),df),
        p_gc=p_gc,q=q,status=reasons)
    diagnostics[[i]] <- data.frame(environment=env_names[i],pc_count=pc,latent_k=cfg$lfmm_k,
        n=nrow(g),residual_df=df,loci=m,tested=sum(is.finite(q)),gc_lambda=inflation,
        significant=sum(q<cfg$alpha,na.rm=TRUE))
}
results <- do.call(rbind,results)
write_tsv(results,"associations.tsv.gz")
write_tsv(results[which(results$q<cfg$alpha),],"discoveries.tsv")
write_tsv(do.call(rbind,diagnostics),"diagnostics.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))
