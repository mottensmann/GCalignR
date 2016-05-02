

ReName = function(Data,ColumnNames=c("ApexRT","StartRT", "EndRT","Area","AreaRel","Height","HeightRel")){
  names(Data) <- ColumnNames
  Data
}


delete_space_colnames <- function(df) {
  names(df) <- str_replace_all(names(df), " ", "")
  df
}

