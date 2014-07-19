using System;
using System.Collections.Generic;

public class SVM {

  private ProblemSetup problem; // Class containing all info
  private int nParticles; 
  private Random r = new Random();

  public SVM(ProblemSetup problem) {
    this.problem = problem;
    nParticles = problem.getNumberOfParticles();
  }

  public matrix generateA() {
    matrix A = new matrix(nParticles - 1, nParticles -1);
    List<double> alphas = new List<double>(); // non-linear parameters
    double min = problem.minGuess(); // Minimum value to look for parameters
    double max = problem.maxGuess(); // Maximum value to look for parameters

    // Prepare randdom values
    for (int i = 0; i < nParticles * (nParticles - 1) / 2; i++) {
      alphas.Add( r.NextDouble() * (max - min) + min);
    }

    // Generate values for every entry in the matrix A
    for (int k = 0; k < nParticles - 1; k++) {
      for (int l = 0; l < nParticles - 1; l++) {
        A[k,l] = A_kl(k,l,alphas);
      }
    }

    return A;
  }

  // Generate the value in A belonging to spot A[k,l]
  private double A_kl(int k, int l, List<double> alphas) {
    double result = 0;

    matrix U_inv = problem.getUInverse(); // Get the inverse U matrix
    List<double> products = new List<double>();
    double B_ijk = 0;
    double B_ijl = 0;
    int count = 0;

    for (int i = 0; i < nParticles; i++) {
      for (int j = 0; j < nParticles; j++) {
        if ( j > i) {
          B_ijk = U_inv[i,k] - U_inv[j,k];
          B_ijl = U_inv[i,l] - U_inv[j,l];
          products.Add(B_ijk * B_ijl);
          count++;
        }
      }
    }

    // Calculate the element based on random values and U_inv matrix
    for (int i = 0; i < count; i++) {
      result += alphas[i] * products[i];
    }

    return result;
  }
}
