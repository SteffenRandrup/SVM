using System;

class CholeskyDecomposition
{
  public matrix L;

  public CholeskyDecomposition(matrix B) {
    matrix A = B.copy();
    /* int size = A.cols; */
    /* L = new matrix (size,size); */
    /* matrix D = new matrix (size,size); */

    /* for (int i = 0; i < size -1; i++) { */
    /*   for (int j = 1; j < size - i; j++) { */
    /*     double ratio = A[i+j,i] / A[i,i]; */
    /*     for (int k = 0; k < size; k++) { */
    /*       A[i+j, k] -= ratio * A[i,k]; */
    /*     } */
    /*     L[i+j,i] = ratio; */
    /*   } */
    /* } */

    /* for (int l = 0; l < size; l++) { */
    /*   L[l,l] = 1; */
    /*   if (A[l,l] < 0) { */
    /*     throw new System.DivideByZeroException(); */
    /*   } */
    /*   D[l,l] = Math.Sqrt(A[l,l]); */
    /* } */
    /* /1* Console.WriteLine("A = "); *1/ */
    /* /1* A.print(); *1/ */
    /* /1* Console.WriteLine("D = "); *1/ */
    /* /1* D.print(); *1/ */
    /* L = L * D; */
  /* } */
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
