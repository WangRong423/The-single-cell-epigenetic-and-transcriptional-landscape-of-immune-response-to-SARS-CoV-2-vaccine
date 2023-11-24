##修改策略：先将zero inflation cutoff值设置得比较小，然后当第一个cluster达到smallest.cluster.num条件，被加入黑名单是，将它对应的zero inflation赋值给zero inflation cutoff
selectCluster_to_proceed_inflation_JXY <- function(inflation.list, IDs, cluster.size.cutoff = 100,blacklist){
  # inflation.list=inflation.tracking
  # IDs=next_round_IDs
  # cluster.size.cutoff = smallest.cluster.num
  blacklist=unique(blacklist)
  cluster.sizes <- table(IDs)
  passed.clusters <- which(cluster.sizes >= cluster.size.cutoff)
  if (!is.null(blacklist)){
    passed.clusters <- passed.clusters[!(passed.clusters%in%blacklist)]
  }
  go_with_higher_inflation <- which.max(inflation.list[passed.clusters])
  selected.cluster <- passed.clusters[go_with_higher_inflation]
  return(selected.cluster)
}




library(JXYlightHippo,lib.loc = "/home/jiangxinyu/R/x86_64-pc-linux-gnu-library/4.1")
library(Seurat,lib.loc = "/home/jiangxinyu/R/x86_64-pc-linux-gnu-library/4.1")
library(sys,lib.loc = "/home/jiangxinyu/R/x86_64-pc-linux-gnu-library/4.1")
library(Signac,lib.loc = "/home/wangrong/R/x86_64-pc-linux-gnu-library/4.2")
time1 <- Sys.time()
#dat <- as_matrix(Gex)#单样本数据用来测试
dat <- readRDS("/database/wangrong/Results/0712_ATAC+RNA/Signac/co_joint/5.14_merge_cojoint.rds")
dat <- dat@assays[["RNA"]]@counts
dat <- as_matrix(dat)
K.round = 600
initial.labels = NULL
initial.round = 0
stop_at = 500
correctByK = FALSE
override.Zscore.cutoff = NULL
smallest.cluster.num = 400
random.num = 5000
move.by.inflation = TRUE
inflation_cutoff <- 0
#lightHIPPO <- function(dat, K.round = 10, initial.labels = NULL, initial.round = 0, stop_at = 500, correctByK = FALSE, override.Zscore.cutoff = NULL, smallest.cluster.num = 200, random.num = 2500, move.by.inflation = TRUE){

  require(irlba)
  total.num.cell <- ncol(dat)
  total.num.gene <- nrow(dat)

  if(!is.null(initial.labels)){
    if(length(initial.labels) != total.num.cell){
      stop("Length of initial group labels doesn't match the number of cell.")
    }
    initial.round <- 0
  }

  if(!is.null(override.Zscore.cutoff)) {
    Zscore.cutoff <- override.Zscore.cutoff
  } else if(correctByK == FALSE){
    Zscore.cutoff <- cut_off_zscore(total.num.gene)#根据基因数计算Zscore阈值 4.54
  } else {
    Zscore.cutoff <- cut_off_zscore(total.num.gene*K)
  }

  if(move.by.inflation == TRUE){

    ### calculate the inflation number for each cluster based on a random set of genes ###
    set.seed(1234567)
    #从1：总gene.num随机出random.num个数
    randomIDs <- sample(1:total.num.gene, random.num)

    if(is.null(initial.labels) & initial.round > 0) {

      initial_clusters <- initialize_HIPPO(dat, initial.round = initial.round, stop_at = stop_at, Zscore.cutoff = Zscore.cutoff)
      next_round_IDs <- initial_clusters$next_round_IDs
      res <- initial_clusters

      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL
      inflation.tracking <- NULL
      for(i in 1:c(initial.round+1)){
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%i], Zscore.cutoff = Zscore.cutoff))
      }
      names(inflation.tracking) <- 1:c(initial.round + 1)

      for(i.round in (initial.round + 1):K.round){

        go_with_higher_inflationID <- selectCluster_to_proceed_inflation_JXY(inflation.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
        selected.dat <- dat[, next_round_IDs%in%go_with_higher_inflationID]
        selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
        selected.ID <- selected.res$selected
        selected.Zscore <- selected.res$Zscore

        new.subset.dat <- selected.dat[selected.ID, ]
        clusterID <- run_kmeans_clustering(new.subset.dat)
        next_round_IDs[next_round_IDs%in%go_with_higher_inflationID][clusterID == 2] <- i.round + 1
        res$sequence <- c(res$sequence, go_with_higher_inflationID)

        inflation.tracking[go_with_higher_inflationID] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%go_with_higher_inflationID], Zscore.cutoff = Zscore.cutoff)
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%(i.round+1)], Zscore.cutoff = Zscore.cutoff))
        names(inflation.tracking)[i.round+1] <- i.round+1

        selected.gene.list[[i.round]] <- selected.ID
        selected.gene.Zscore[[i.round]] <- selected.Zscore
      }
      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"
      res$initial.clusters <- NULL

    } else if(is.null(initial.labels) & initial.round == 0) {

      res <- NULL
      res$sequence <- NULL
      inflation.tracking <- NULL
      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL
      blacklist <- NULL
      i.label <- 0
      for(i.round in (initial.round + 1):K.round){
        i.round <- i.round-i.label
        if(i.round == 1){

          selected.res <- select_features_full(dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          selected.Zscore <- selected.res$Zscore
          new.subset.dat <- dat[selected.ID, ]
          clusterID <- run_kmeans_clustering(new.subset.dat)
          next_round_IDs <- clusterID
          #随机抽样检查每一簇中的inflation gene
          inflation.tracking[1] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs==1], Zscore.cutoff = Zscore.cutoff)
          inflation.tracking[2] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs==2], Zscore.cutoff = Zscore.cutoff)
          names(inflation.tracking) <- c(1:2)

          selected.gene.list[[i.round]] <- selected.ID
          selected.gene.Zscore[[i.round]] <- selected.Zscore

        } else {
          #选出inflation cluster ID
          print(inflation.tracking)
          go_with_higher_inflationID <- selectCluster_to_proceed_inflation_JXY(inflation.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num,blacklist = blacklist)
          if(!length(go_with_higher_inflationID)){
            break
          }
          selected.dat <- dat[, next_round_IDs%in% go_with_higher_inflationID]
          selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          selected.Zscore <- selected.res$Zscore

          new.subset.dat <- selected.dat[selected.ID, ]
          #将选择出的inflation cluster 继续分
          clusterID <- run_kmeans_clustering(new.subset.dat)
          #table(clusterID)
          if(length(blacklist)==1){
            inflation_cutoff <- inflation.tracking[blacklist]
          }
          blacklist <- c(blacklist,which(inflation.tracking<0.8*inflation_cutoff))
          blacklist <- unique(blacklist)
          flag=table(clusterID)<smallest.cluster.num
          if (flag[1]|flag[2]){
            #加入黑名单 表明这个cluster无需再分
            blacklist <- c(blacklist,go_with_higher_inflationID)
            i.label <- i.label+1
            next
          }
          #将next_round_IDs中的inflation cluster分出的“2”贴上新的标签
          next_round_IDs[next_round_IDs%in%go_with_higher_inflationID][clusterID == 2] <- i.round + 1
          print(table(next_round_IDs))
          #将inflationID做一个记录
          res$sequence <- c(res$sequence, go_with_higher_inflationID)
          #更新inflation.tracking（原来inflation cluster的inflation number会下降）
          inflation.tracking[go_with_higher_inflationID] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%go_with_higher_inflationID], Zscore.cutoff = Zscore.cutoff)
          #check new cluster的inflaton gene number  and lables
          inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%(i.round+1)], Zscore.cutoff = Zscore.cutoff))
          names(inflation.tracking)[i.round+1] <- i.round+1

          selected.gene.list[[i.round]] <- selected.ID
          selected.gene.Zscore[[i.round]] <- selected.Zscore
          print(i.round)
        }

      }
      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"
      res$initial.clusters <- NULL

    } else if(!is.null(initial.labels)) {

      res <- NULL
      res$sequence <- NULL
      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL

      next_round_IDs <- initial.labels

      inflation.tracking <- NULL
      for(i in 1:max(initial.labels)){
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%i], Zscore.cutoff = Zscore.cutoff))
      }
      names(inflation.tracking) <- 1:max(initial.labels)

      for(i.round in (max(initial.labels)):(K.round+max(initial.labels)-1)){

        go_with_higher_inflationID <- selectCluster_to_proceed_inflation(inflation.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
        selected.dat <- dat[, next_round_IDs%in% go_with_higher_inflationID]
        selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
        selected.ID <- selected.res$selected
        selected.Zscore <- selected.res$Zscore

        new.subset.dat <- selected.dat[selected.ID, ]
        clusterID <- run_kmeans_clustering(new.subset.dat)
        next_round_IDs[next_round_IDs%in%go_with_higher_inflationID][clusterID == 2] <- i.round + 1
        res$sequence <- c(res$sequence, go_with_higher_inflationID)

        inflation.tracking[go_with_higher_inflationID] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%go_with_higher_inflationID], Zscore.cutoff = Zscore.cutoff)
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%(i.round+1)], Zscore.cutoff = Zscore.cutoff))
        names(inflation.tracking)[i.round+1] <- i.round+1

        selected.gene.list[[i.round]] <- selected.ID
        selected.gene.Zscore[[i.round]] <- selected.Zscore

      }

      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+max(initial.labels)
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Truncated"
      res$initial.clusters <- initial.labels
    }

  }

  if(move.by.inflation == FALSE) {

    ### calculating first 10 PCs
    first_pcs <- tryCatch(expr = {
      irlba::irlba(log1p(dat), 10)$v
    }, error = function(e) NA, warning = function(w) NA)

    if(is.null(initial.labels) & initial.round > 0) {

      initial_clusters <- initialize_HIPPO(dat, initial.round = initial.round, stop_at = stop_at, Zscore.cutoff = Zscore.cutoff)
      next_round_IDs <- initial_clusters$next_round_IDs
      res <- initial_clusters

      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL
      cluster.heterogeneity.tracking <- NULL

      kkk <- apply(first_pcs, 2, function(x){
        tapply(x, next_round_IDs, var)
      })
      cluster.heterogeneity.tracking <- apply(kkk, 1, sum)*table(next_round_IDs)

      for(i.round in (initial.round + 1):K.round){

        go_with_larger_varianceID <- selectCluster_to_proceed(cluster.heterogeneity.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
        selected.dat <- dat[, next_round_IDs%in%go_with_larger_varianceID]
        selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
        selected.ID <- selected.res$selected
        selected.Zscore <- selected.res$Zscore

        new.subset.dat <- selected.dat[selected.ID, ]
        clusterID <- run_kmeans_clustering(new.subset.dat)

        new.clusters.var <- apply(first_pcs[next_round_IDs%in%go_with_larger_varianceID, ], 2, function(x){
          tapply(x, clusterID, var)
        })

        next_round_IDs[next_round_IDs%in%go_with_larger_varianceID][clusterID == 2] <- i.round + 1
        res$sequence <- c(res$sequence, go_with_larger_varianceID)

        new.cluster.heterogeneity.tracking <- apply(new.clusters.var, 1, sum)*table(clusterID)
        cluster.heterogeneity.tracking[go_with_larger_varianceID] <- new.cluster.heterogeneity.tracking[1]
        cluster.heterogeneity.tracking <- c(cluster.heterogeneity.tracking, new.cluster.heterogeneity.tracking[2])
        names(cluster.heterogeneity.tracking)[i.round + 1] <-  i.round + 1

        selected.gene.list[[i.round]] <- selected.ID
        selected.gene.Zscore[[i.round]] <- selected.Zscore

      }

      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"
      res$initial.clusters <- NULL

    }


  }
  time2 <- Sys.time()
  print(time2-time1)
  saveRDS(res,"/database/jiangxinyu/Result/Data/0516_Total_cojoint_res_400cells_random5000_improve0.8.rds")
  ##plot
  # summarizing_dat <- summarize_current_zero_proportions(dat[randomIDs, ], res$next_round_IDs)
  # plot_dat_per_cluster_inflation <- visualize_current_zero_proportions(summarizing_dat)
  #
  # library(HIPPO)
  # Gex.sce = SingleCellExperiment(assays = list(counts = Gex))
  # hippo_diagnostic_plot(Gex.sce,
  #                       show_outliers = T,
  #                       zvalue_thresh = 20)

# return(res)

#}
