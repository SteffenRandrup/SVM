using System;

class CholeskyDecomposition
{
  public matrix L;

  public CholeskyDecomposition(matrix A) {
    int size = A.cols;
    L = new matrix (size,size);
    matrix D = new matrix (size,size);

    for (int i = 0; i < size -1; i++) {
      for (int j = 1; j < size - i; j++) {
        double ratio = A[i+j,i] / A[i,i];
        for (int k = 0; k < size; k++) {
          A[i+j, k] -= ratio * A[i,k];
        }
        L[i+j,i] = ratio;
      }
    }

    for (int l = 0; l < size; l++) {
      L[l,l] = 1;
      D[l,l] = Math.Sqrt(A[l,l]);
    }
    L = L * D;
  }

}
