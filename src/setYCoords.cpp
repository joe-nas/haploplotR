#include <Rcpp.h>
using namespace Rcpp;

//' internal functions used in genesGrob to specify y position of transcripts
//' @name setYCoordsFUN
//' @param txs GRanges or IRanges object specifying start and end postitions of transcripts
//' @return integer vector indicating y postions of transcripts
//' @examples
//' where <- GRanges(Rle('chr1'), IRanges(start = 209789729, end = 210189729))
//' COMPLETE EXAMPLE

// [[Rcpp::export]]
DataFrame setYCoordsFUN(DataFrame df){

  // access the columns
  IntegerVector start = df["start"];
  IntegerVector end = df["end"];

  int nrow = start.length();

  // new column y
  IntegerVector y(nrow);

  //algo
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < i; j++){
      if(start[i] <= end[j] && start[j] <= end[i]){
        y[i] = y[i] + 1;
      }
    }
  }

  for(int i = 0; i < nrow; i++){
    y[i] = y[i] + 1;
  }

  // return a new data frame
  return DataFrame::create(df, _["y"]=y);
}
