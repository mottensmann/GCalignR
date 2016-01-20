RetentionCutoff <- function(Data, Low=8, High=NULL){
  # RetentionCutoff removes all Retention Times below the Threshold specified by Low (default 8s).
  # In addition Retention Times above a time defined by the Value of High (Default is Null)
  # can be applied.
    if (is.null(High)){
        out <- subset(Data, ApexRT > Low)
    }else{
        out <- subset(Data, ApexRT > Low & ApexRT< High)
    }
    out
}
