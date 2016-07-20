function_call <- function(call,FUN="align_chromatograms"){
form <- formals(FUN)
for ( n in names(form)){ # for every args of the function
    if (!(n %in% names(call))){ # Find args not called
        call <- append(call,form[n])  ## add missing args
    }
}
type <- as.vector(which(lapply(call, function(x) out <- class(x))!="NULL"))
call[type] <- lapply(call[type], function(x) x <- as.character(x))
call[-type] <- lapply(call[-type], function(x) x <- "NULL")
call <- do.call(rbind,call)
call <- t(as.data.frame(call))
row.names(call) <- NULL
return(call)
}

