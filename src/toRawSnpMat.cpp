#include <Rcpp.h>
using namespace Rcpp;

//' Create matrix with type Raw from matrices describing heterogeneity and presence of the alternative allele
//' @name toRawSnpMat
//' @param ishet Integer matrix describing heterogeneity
//' @param hasalt Integer matrix describing whether or not the alternaive allele is present.
//' @return matrix A matrix of encoding the genotype as type Raw. 01 homozygous reference allele, 02 heterozygous, 03 homozygous alternative allele, 04 missing and else.

// [[Rcpp::export]]
RawMatrix toRawSnpMat(IntegerMatrix ishet, IntegerMatrix hasalt){
  int nrow = hasalt.nrow();
  int ncol = hasalt.ncol();

  int a;
  int h;

  RawMatrix out(nrow, ncol);

  for(int i = 0; i < nrow; i++ ){
    for(int j = 0; j < ncol; j++ ){
      a = hasalt(i,j);
      h = ishet(i,j);

      if(a == 0 && h == 0){
        out(i,j) = 1;
      }else if(h == 1){
        out(i,j) = 2;
      }else if(h == 0 && a == 1){
        out(i,j) = 3;
      }else{
        out(i,j) = 4;
      }
    }
  }

  return out;
}
