using System;

class CholeskyDecomposition
{
  public matrix L;

  public CholeskyDecomposition(matrix B) {
    matrix A = B.copy();

    int size = A.cols;
    L = new matrix(size, size);
    for (int i = 0; i < size; i++) {
      for (int k = 0; k < (i+1); k++) {
        double sum = 0;
        for (int j = 0; j < k; j++) {
          sum += L[i,j] * L[k,j];
        }
        if ( i == k) {
          L[i,k] = Math.Sqrt(A[i,i] - sum);
        } else {
          L[i,k] = 1 / L[k,k] * (A[i,k] - sum);
        }
      }
    }
}
}
