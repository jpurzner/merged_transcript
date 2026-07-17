plot_group_timeline <- function (expr_list, meta_list, gene2cluster, p_time = NULL, 
                                     ylimits = c(-2,2), include_n = TRUE, col_num = 4, 
                                     label_len = 20, text_size = 16, order_clust = FALSE, 
                                     ribbon = TRUE, x_label = "none",  add_t50 = NULL, time_ref = NULL) {
  
  #preseves the order of the factors unless order_clust = TRUE
  # x_label NULL, "all", "drop" , "remove"
  
  # t50_to_compressed() is a sibling helper in R/lib. It is normally sourced
  # by mt_source_lib(); source it from this file's directory as a fallback so
  # plot_group_timeline() is self-contained wherever the repo lives.
  if (!exists("t50_to_compressed", mode = "function")) {
    .self_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
    if (!is.null(.self_dir) && file.exists(file.path(.self_dir, "t50_to_compressed.R")))
      source(file.path(.self_dir, "t50_to_compressed.R"))
  }
  require(ggplot2)
  require(plyr)
  require(ggrepel)
  require(reshape2)
  require(dplyr)
  
  #require(naturalsort)
  
  # loop through the list of expression objects and extract data, merge with community clusters, 
  #  flatten the data and rbind into a common df for use with ggplot2
 
#print
  original_gene2cluster <- gene2cluster
  
  toDE <- function(DEin) {
    DE <- as.data.frame(DEin, row.names = names(DEin))
    DE$genes <- row.names(DE)
    row.names(DE) <- c(1:(dim(DE)[1]))
    DE <- DE[,c(2,1)]
    colnames(DE)[2] <- "cluster_id"
    colnames(DE)[1] <- "gene"
    return(DE)
  }  
  
  # handle the case where the row.names are gene names
  if(is.null(dim(gene2cluster)[2])) {   
      gene2cluster <- toDE(gene2cluster)
      #print(head(gene2cluster))
  }
  else if (dim(gene2cluster)[2] == 1 )  {
    gene2cluster <- toDE(gene2cluster)
    #print(head(gene2cluster))
  }
  
  if (is.factor(gene2cluster)) {
    gene2cluster$cluster_id <- factor(gene2cluster$cluster_id, levels = levels(gene2cluster))
    
  }
  
  #print(head(gene2cluster))  
  colnames(expr_list[[1]]) <- c(1:(dim(expr_list[[1]])[2]))  
  expr <- merge(x = gene2cluster, y = expr_list[[1]], by.x = 1, by.y = 0)
  
  expr <- melt(expr, id.vars = c("gene", "cluster_id"))
  #print(expr)
  colnames(expr)[3] <- "time_idx"
  expr$ds <- names(meta_list)[1]
  expr$time <- sapply(expr$time_idx, function(x) meta_list[[1]][x])
  expr <- expr[,c(1,2,6,3,5,4)]
  #print(head(expr))
  #print(tail(expr))
  #print(length(expr_list))

  
  
  if (length(expr_list) > 1) {
    for (i in 2:length(expr_list)) {
      expr_tmp <- expr_list[[i]]
      expr_tmp <- merge(x = gene2cluster, y = expr_tmp, by.x = 1, by.y = 0)
      expr_tmp <- melt(expr_tmp, id.vars = c("gene","cluster_id"))
      colnames(expr_tmp)[3] <- "time_idx"
  
      expr_tmp$ds <- names(meta_list)[i]
      expr_tmp$time <- sapply(expr_tmp$time_idx, function(x) meta_list[[i]][x])
      expr_tmp <- expr_tmp[,c(1,2,6,3,5,4)]
    
      expr <- rbind(expr, expr_tmp)
    } # end for loop 
  } # end if 
  
  expr$time <- factor(expr$time)

  # to overlap plot the data we need to use a pseudotime 
  # this can be specified using the following df 
  # 
  # ds  time   pseudotime
  #Scott_GNP  E15  0
  #Scott_GNP  P1  0 
  
  # series of re-ordering to make the pseudotime df 
  ### this block of code is defunct unless we need to seperate different ds with same time
  #pseudo_time <- ddply(expr, c("ds", "time"), summarise, N = length(value))
  #pseudo_time$N <- NULL
  #pseudo_time$ds <- factor(pseudo_time$ds, levels = ds_names)
  #pseudo_time <- pseudo_time[order(pseudo_time$ds),]
  #pseudo_time$time <- factor(pseudo_time$time, levels = p_time)
  #pseudo_time <- pseudo_time[order(pseudo_time$time),]
  #pseudo_time$ptime <-  as.numeric(pseudo_time$time)  
  #row.names(pseudo_time) <- seq(length=nrow(pseudo_time)) 
  
  ## TODO automatically try and make the pseudotime table
  # finds the data set with the highest number of time points
  #ds_largest <- table(pseudo_time$ds)
  #ds_largest <- names(ds_largest[order(ds_largest, decreasing  = T)])[1]
  #largest_time <- pseudo_time[pseudo_time$ds == ds_largest, 2]
  
  ##END of defunct code 
  
  # much shorter pseudo time setting 
    
  
  # if pseudotime provided 
  if (!is.null(p_time)) {
    # check that pseudo time names match 
    ptime_list <- seq(length=length(p_time))
    names(ptime_list) <- p_time
    ptime_list <- factor(ptime_list)
    ptime_list <- as.data.frame(ptime_list)  
  } else {
    if (length(meta_list) > 1) {
      ptime_list = as.data.frame(unique(as.vector(unlist(lapply(meta_list, levels)))))  
    } else {
      ptime_list = as.data.frame(levels(meta_list[[1]]))
      #print("test")
    } 
  }

  #print(head(expr))
  #print(ptime_list)
  expr <- merge(x = expr, y = ptime_list, by.x = 3, by.y = 0, all.x = T)
  
  expr <- expr[,c(2:5,1,7,6)]
  colnames(expr)[6] <- "ptime"
  expr$ptime <- factor(expr$ptime)
  expr$cluster_id <- factor(expr$cluster_id)
  # add the pseudo to df
  #expr$ptime <- ptime_list[expr$time]
  #expr <- expr[,c(1:5,7,6)]
  #print(head(expr))
  
  #print(head(expr))

  
  # collapse the data frame to get averages
  mean_expr <- ddply(expr, c("cluster_id", "ds", "time", "ptime"), summarise, mean = mean(value), sd = sd(value), sem = sd(value)/sqrt(length(value)))
  #mean_expr$ptime <- factor(mean_expr$ptime)
  
  if (order_clust) {
    new_level <- levels(mean_expr$cluster_id)[suppressWarnings(order(as.numeric(levels(mean_expr$cluster_id))))]
    mean_expr$cluster_id <- factor(mean_expr$cluster_id, levels = new_level)
    #print("number")
  }
  
  #print(head(mean_expr))
  #print(levels(mean_expr$cluster_id))
  
  clust_num <- as.data.frame(table(gene2cluster$cluster_id))
  colnames(clust_num) <- c("cluster_id", "clust_num")
  
re_name <- {}
  if (include_n) {
    clust_num$title <- paste(clust_num$cluster_id,"n =",clust_num$clust_num)   
    mean_expr <- merge(x = mean_expr, y = clust_num, by.x = 1, by.y = 1)
    mean_expr$title <- factor(mean_expr$title)
    re_name <- mean_expr[!duplicated(mean_expr$cluster_id),c("cluster_id", "title")]
    level_vals <- re_name$title[match(levels(mean_expr$cluster_id), re_name$cluster_id)]
    mean_expr$cluster_id <- factor(mean_expr$title, levels = level_vals)
  } else { 
    # maintain factor order 
    first_ind <- which(!duplicated(mean_expr$cluster_id))
    names(first_ind) <- mean_expr$cluster_id[first_ind]
    first_ind <- first_ind[match(names(first_ind), levels( mean_expr$cluster_id))]
    mean_expr$cluster_id <- factor(mean_expr$cluster_id, levels = names(first_ind))
  }
  
  
  old_levels = factor(gsub("_", " ", levels(mean_expr$cluster_id)))   
  mean_expr$cluster_id <- factor(gsub("_", " ", mean_expr$cluster_id), levels = old_levels)
  
  # merge t50 by cluster_id 
  if (!is.null(add_t50)) {
    add_t50 <- merge(x = add_t50,  y = gene2cluster, by.x = 1, by.y = 1, all = FALSE)
    add_t50 <- add_t50 %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(median_avg_t50 = median(avg_t50))
    # translate the 
    add_t50$median_avg_t50_cor <-   sapply(add_t50$median_avg_t50, function (x) t50_to_compressed(x, time_ref ))
    
    # modify cluster_id with include_n 
    if (include_n) {
      add_t50$cluster_id <- as.character(add_t50$cluster_id)
      add_t50$cluster_id <- plyr::mapvalues(add_t50$cluster_id,  from = as.character(re_name$cluster_id), to = as.character(re_name$title))
      add_t50$cluster_id 
    }
    
    add_t50$cluster_id <- factor(gsub("_", " ", add_t50$cluster_id), levels = levels( mean_expr$cluster_id))
    
   }
  

  
  print(head(add_t50))
 
  print(head(mean_expr))
  #print(head(mean_expr))
  #print(levels(mean_expr$cluster_id))
  
  #clust_name_dict <-clust_num$title
  #names(clust_name_dict)  <- clust_num$cluster_id
  
  
  #print(head(mean_expr))
  #print(table(mean_expr$cluster_id))
  
  p <- ggplot(mean_expr, aes(x=as.numeric(factor(ptime)), y=mean, colour = ds, group = ds)) + 
    geom_line(size = 1.25) 
  
  if (ribbon) {
    p <- p + geom_ribbon(aes(ymax = mean + sd, ymin = mean - sd, fill = ds), alpha = 0.3, colour=NA) 
  }
  p <- p + geom_point() +
    facet_wrap(~cluster_id, labeller=labeller(cluster_id = label_wrap_gen(label_len)), ncol = col_num, scales ="free") 
  

  
  # handle case with a list of genes      
  if (length(expr_list) > 1) {
    p <- p + geom_text_repel(aes(label=time), 
                    nudge_y = ifelse(mean_expr$ds == names(meta_list)[2], -0.5, ifelse(mean_expr$ds == names(meta_list)[1], 0.5, 0.75)), 
                    size = text_size*(4/14), force =1 ) 
  } else {
    p <- p + geom_text_repel(aes(label=time), 
                      nudge_y = 0.5, 
                      size = text_size*(4/14), force =1 ) 
  }  
  
  # add t50 line
  if (!is.null(add_t50)) {
    p <- p + geom_point(data = add_t50, aes(x = median_avg_t50_cor, y = 0), 
                        colour = "black", 
                        fill = "black",  
                        inherit.aes = FALSE ) + 
      geom_segment(data = add_t50, aes(x= median_avg_t50_cor, xend  = median_avg_t50_cor, y = 0, yend = -Inf ),
                   inherit.aes = FALSE) + 
      
      geom_text(data = add_t50, aes( x = median_avg_t50_cor + 0.75 , y = ylimits[1] + 0.25, label = "t50"), 
                inherit.aes = FALSE, 
                size = text_size * 5/14 )  
    
    
  }


  p <- p + coord_cartesian(ylim=ylimits) +
    ylab("log2 fold change") 

  
  p <- p + theme(legend.position="bottom", 
                 axis.line = element_line(size = 1),
                 legend.title=element_blank(), 
                 panel.border=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 plot.background=element_blank(),
                 panel.background=element_blank(),
                 strip.background = element_blank(),
                 text=element_text(size=text_size), 
                 strip.text = element_text(size=text_size))
  
  if (x_label == "all") {
    p <- p + xlab("developmental time") + 
      scale_x_continuous(labels = p_time,  breaks = c(1:length(p_time))) + 
      theme(axis.text.x = element_text(size=text_size-4)) 
  } else if (x_label == "drop") {
    # stagger the x-axis labels by adding a newline alternatively 
    # from https://stackoverflow.com/questions/19567986/overlapping-axis-labels-in-r
    alt_axis <- function( labels ) { 
      fixedLabels <- c() 
      for ( l in 1:length( labels ) ) { 
        fixedLabels <- c( fixedLabels, paste0( ifelse( l %% 2 == 0, '', '\n' ), labels[l] ) ) 
        } 
      return( fixedLabels ) 
    } 
    
    p <- p + xlab("developmental time") + 
      scale_x_continuous(labels = alt_axis(p_time), breaks = c(1:length(p_time))) + 
      theme(axis.text.x = element_text(size=text_size-4)) 
  } else if (x_label == "remove") {
    # stagger the x-axis labels by adding a newline alternatively 
    # from https://stackoverflow.com/questions/19567986/overlapping-axis-labels-in-r
    alt_axis <- function( labels ) { 
      fixedLabels <- c() 
      for ( l in 1:length( labels ) ) { 
        fixedLabels <- c( fixedLabels, ifelse( l %% 2 == 0, '',  labels[l] ) ) 
      } 
      return( fixedLabels ) 
    } 
    
    p <- p + xlab("developmental time") + 
      scale_x_continuous(labels = alt_axis(p_time), breaks = c(1:length(p_time))) + 
      theme(axis.text.x = element_text(size=text_size-4)) 
  } else {
    p <- p + theme(axis.text.x =element_blank(),
                   axis.title.x=element_blank(),
                   axis.ticks.x =element_blank()) 
  }
    
  p
  
  return(p)
  
  # link genes to co_id
  
  
  
  
  
} 