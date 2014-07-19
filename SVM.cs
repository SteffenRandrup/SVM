using System;
using System.Collections.Generic;

public class SVM {

  private ProblemSetup problem; // Class containing all info
  private int nParticles;
  private Random r = new Random();

  public SVM(ProblemSetup problem) {
    this.problem = problem;
    nParticles = problem.getNumberOfParticles();
    generateB().print();
  }

  // Generate matrix containing non-linear parameters for the gaussian test
  // functions transformed to center of mass system
  private matrix generateA() {
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

  // Symmetrize a matrix so that there is anti symmetry on swapping
  // indistinguishable fermions and symmetry on swapping indistinguishable
  // bosons
  private matrix symmetrize(matrix A) {
    // Generate permutations of particles
     List<int> perm = new List<int>();
     for (int i = 0; i < nParticles; i++){
       perm.Add(i);
     }

     // Make a list with 1 permutation. It will be used for getting permutations
     // of
     List<Permutation> ll = new List<Permutation>();
     ll.Add(new Permutation(perm,1,-1)); // fermions

     matrix result = new matrix(A.rows,A.cols);

     foreach(Permutation ps in perms(ll)) {
       // Generate C_i according to the permutation and transform to jacobi
       matrix Ci = generateC(ps.ToArray());
       matrix tmp = problem.getUInverse() * Ci * problem.getU();
       matrix Pi = new matrix(A.rows, A.cols);

       // Remove last row and column from matrix
       for(int i = 0; i < A.rows; i++) {
         for(int j = 0; j < A.cols; j++) {
          Pi[i,j] = tmp[i,j];
         }
       }

       // There might be an issue with size here
       result += ps.Parity() * (Pi.transpose() *A * Pi);
     }

     return 1/Math.Sqrt(Misc.factorial(A.rows)) * result;
  }

  // Recursively generate a list of permutations.
  // The permutation class contains the order and parity of a
  // permutation.
  private List<Permutation> perms(List<Permutation> permlist) {
    List<Permutation> result = new List<Permutation>();
    foreach (Permutation p in permlist) {
      // If length is 1 there is only one permutation
      if ( p.Length() == 1) {
        result.Add(p.Clone());}
      // If length is 2 swap the content
      else if ( p.Length() == 2) {
        result.Add(p.Clone());
        result.Add(p.swap(0,1));}
      // For every number in permutation remove it and generate
      // permutation of the remaining. Then append them to the number
      else {
        for (int i = 0; i < p.Length(); i++) {
          Permutation pc = p.Clone();
          Permutation pp = pc.Pop(i);
          List<Permutation> list = new List<Permutation>();
          list.Add(pc);
          foreach(Permutation put in perms(list)) {
            result.Add(pp.Clone().Append(put));
          }
        }
      }
    }
    return result;
  }

  // Generate a matrix representing the given permutation
  private matrix generateC(int[] permutation) {
    matrix C = new matrix(nParticles,nParticles);
    for(int i = 0; i < nParticles; i++) {
      for(int j = 0; j < nParticles; j++) {
        if(j == permutation[i]) {
          C[i,j] = 1;
        }
      }
    }
    return C;
  }

  // Calculate matrix element for H matrix
  private double matrixElement(matrix A, matrix B) {
    return kineticEnergyElement(A,B) + potentialEnergyElement(A,B);
  }

  // Calculate the overlap of two test functions represented by matrices
  private double overlapElement(matrix A, matrix B) {
    QRdecomposition qr = new QRdecomposition(A+B);
    double determinant = 1;
    for (int i = 0; i < A.cols; i++) {
      determinant *= qr.R[i,i];
    }
    double nominator = Math.Pow( 2 * Math.PI, A.cols);
    return Math.Pow(nominator / determinant, 3/2);
  }

  // Calculates the matrix element for the kinetic energy
  private double kineticEnergyElement(matrix A, matrix B) {
    QRdecomposition qr = new QRdecomposition(A+B);
    matrix inverse = qr.inverse();
    return 3 * (A * inverse * B * problem.getLambda()).trace() * overlapElement(A,B);
  }

  // Calculate the matrix element of the potential energy
  // Currently only coulomb potential is supported
  private double potentialEnergyElement(matrix A, matrix B) {
    return coulombPotentialEnergy(A,B);
  }

  // Calculate matrix element of a coulomb potential given the charges
  // of the elements
  private double coulombPotentialEnergy(matrix A, matrix B) {
    double epsilon_0 = 1; // Unit is set to 1 maight want to make static units
    int N = problem.getNumberOfParticles() - 1; // For easier readability
    double result = 0;


    matrix inverse = new QRdecomposition(A+B).inverse();
    // For every particle pair (only counting once) calculate interaction
    List<Particle> particles = problem.getParticles();
    for (int i = 1; i < N; i++) {
      for(int j = 0; j < i; j++) {
        double p_ij = 0;
        matrix Uinv = problem.getUInverse();
        for (int k = 0; k < N; k++) {
          for (int l = 0; l < N; l++) {
            double B_ijk = Uinv[i,k] - Uinv[j,k];
            double B_ijl = Uinv[i,l] - Uinv[j,l];
            p_ij += B_ijk * inverse[k,l] * B_ijl;
          }
        }
        Console.WriteLine("Pijk = " + p_ij);

        result += particles[i].getCharge() * particles[j].getCharge() / (4 * Math.PI * epsilon_0) * Math.Pow( 1 / (2 * p_ij * Math.PI), 3/2) * overlapElement(A,B);
      }
    }

    return result;
  }

  private matrix generateB() {
  
  }
}
