
# extend peak function
extension <- function(peak_df,bp) {
  inds <- peak_df$strand == "+"
  peak_df$end[inds] = peak_df$start[inds] + bp
  peak_df$start[!inds] = peak_df$end[!inds] - bp
  return(peak_df)
}

# transformation function
transformation <- function(peak_df,bp) {
  inds <- peak_df$strand == '+'
  # transform the + strand -50 and +50 from the start site will be the new peak
  peak_df$start[inds] = peak_df$start[inds] - bp
  peak_df$end[inds] = peak_df$start[inds] + bp
  # modify the - strand -50 and +50 from the end site will be the new peak
  peak_df$start[!inds] = peak_df$end[!inds] - bp
  peak_df$end[!inds] = peak_df$end[!inds] + bp
  return(peak_df)
}


# another transformation

extend_grange <- function(pk) {
  tmp_start = ifelse(strand(pk)=='+',start(pk)-101,end(pk)-100)
  tmp_end   = ifelse(strand(pk)=='-',end(pk)+100,start(pk)+99)
  start(pk) = tmp_start
  end(pk)   = tmp_end
  return(pk)
}


extend_grange <- function(pk,bp) {
  tmp_start = ifelse(strand(pk)=='+',start(pk)-(bp+1),end(pk)-(bp))
  tmp_end   = ifelse(strand(pk)=='-',end(pk)+(bp),start(pk)+(bp-1))
  start(pk) = tmp_start
  end(pk)   = tmp_end
  return(pk)
}