time_expr_avg <- function(expr_align_list) {

gene_num <-  dim(expr_align_list[[1]])[1]
expr_avg <- t(sapply(c(1:gene_num), 
                     function(j) apply( sapply(1:length(expr_align_list), 
                                               function(k) expr_align_list[[k]][j,]), 1, 
                                        function (i) mean(i, na.rm = TRUE ))))

row.names(expr_avg) <- row.names(expr_align_list[[1]]) 
return(expr_avg)


}