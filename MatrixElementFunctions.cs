using System;
using System.Collections.Generic;
using math = System.Math;

public class MEF {
  ProblemSetup problem;

  public MEF(ProblemSetup problem) {
    this.problem = problem;
  }

  // Calculate matrix element for H matrix
  public double matrixElement(matrix A, matrix B) {
    /* return kineticEnergyElement(A,B) + potentialEnergyElement(A,B); */
    double kin = kineticEnergyElement(A,B);
    double pot = potentialEnergyElement(A,B);
    /* Console.WriteLine("Kinetic energy element " + kin); */
    /* Console.WriteLine("Potential energy element " + pot); */
    return kin + pot;
  }

  // Calculate the overlap of two test functions represented by matrices
  public double overlapElement(matrix A, matrix B) {
    QRdecomposition qr = new QRdecomposition(A+B);
    double determinant = 1;
    for (int i = 0; i < A.cols; i++) {
      determinant *= qr.R[i,i];
    }
    double nominator = Math.Pow( 2 * Math.PI, A.cols);
    return Math.Pow(nominator / determinant, 3.0/2.0);
  }

  // Calculates the matrix element for the kinetic energy
  public double kineticEnergyElement(matrix A, matrix B) {
    QRdecomposition qr = new QRdecomposition(A+B);
    matrix inverse = qr.inverse();

    matrix tmp = A * inverse * B * problem.getLambda();
    return 3.0/2.0 * tmp.trace() * overlapElement(A,B);
  }

  // Calculate the matrix element of the potential energy
  // Currently only coulomb potential is supported
  public double potentialEnergyElement(matrix A, matrix B) {
    return coulombPotentialEnergy(A,B);
  }

  // Calculate matrix element of a coulomb potential given the charges
  // of the elements
  public double coulombPotentialEnergy(matrix A, matrix B) {
    double epsilon_0 = 1; // Unit is set to 1 maight want to make static units
    /* int N = problem.getNumberOfParticles() - 1; // For easier readability */
    int N = problem.getNumberOfParticles();
    double result = 0;


    matrix inverse = new QRdecomposition(A+B).inverse();
    // For every particle pair (only counting once) calculate interaction
    List<Particle> particles = problem.getParticles();
    matrix Uinv = problem.getUInverse();
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < j; i++) {
        double p_ij = 0;
        for (int k = 0; k < N - 1; k++) {
          for (int l = 0; l < N - 1; l++) {
            double b_ijk = Uinv[i,k] - Uinv[j,k];
            double b_ijl = Uinv[i,l] - Uinv[j,l];
            p_ij += b_ijk * inverse[k,l] * b_ijl;
          }
        }
        /* result +=  particles[i].getCharge() * particles[j].getCharge() / (4 * Math.PI * epsilon_0) * Math.Sqrt(2 / (p_ij * Math.PI)) * overlapElement(A,B); */
        result +=  particles[i].getCharge() * particles[j].getCharge() * 2.0 / epsilon_0 * Math.Sqrt(1/(p_ij * Math.PI * 2)) * overlapElement(A,B);
      }
    }

    return result;
  }
}
