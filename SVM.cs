using System;
using System.Collections.Generic;

public class SVM {

  private ProblemSetup problem; // Class containing all info
  private int nParticles;
  private Random r = new Random();
  private MEF mef;

  public SVM(ProblemSetup problem) {
    this.problem = problem;
    nParticles = problem.getNumberOfParticles();
    mef = new MEF(problem);
  }

  public vector run(int size, int testsize) {
    double E0 = double.PositiveInfinity;
    List<matrix> basis = new List<matrix>();
    vector result = new vector(size);
    for(int i = 0; i < size; i++) {
      List<matrix> testfunctions = generateTestFunctions(testsize);
      List<matrix> bestBasisSoFar = new List<matrix>(basis);;
      for (int j = 0; j < testfunctions.Count; j++) {
        List<matrix> tmpbasis = new List<matrix>(basis);
        tmpbasis.Add(testfunctions[j]);
        /* if (validateB(generateB(tmpbasis))){ */
          vector v = eigenValues(tmpbasis);
          if (v[0] < E0) {
            bestBasisSoFar = tmpbasis;
            E0 = v[0];
          }
        /* } */
      }
      result[i] = E0;
      basis = bestBasisSoFar;
    }
    return result;
  }

  public matrix permute(matrix M, Permutation p) {
    // TODO implement permutating a single matrix according to permu
    return null;
  }


  /* public vector run(int reps) { */
  /*   List<matrix> basis = generateBasis(reps); */
  /*   int count = 1; */
  /*   bool changed = false; */
  /*   double lowestEigen = eigenValues(basis)[0]; */
  /*   List<double> eigens = new List<double>(); */
  /*   eigens.Add(lowestEigen); */
  /*   Console.WriteLine(); */
  /*   Console.WriteLine("Lowest eigen value is " + lowestEigen); */
  /*   Console.WriteLine(); */
  /*   for (int i = 0; i < basis.Count; i++) { */
  /*     while (!changed) { */
  /*       List<matrix> tmpbasis = new List<matrix>(basis); */
  /*       tmpbasis[i] = generateTestFunctions(1)[0]; */
  /*       if (validateB(generateB(tmpbasis))) { */
  /*         double energy = eigenValues(tmpbasis)[0]; */
  /*         if (energy < lowestEigen) { */
  /*           changed = true; */
  /*           lowestEigen = energy; */
  /*           basis = tmpbasis; */
  /*           Console.WriteLine("Changed"); */
  /*           eigens.Add(energy); */
  /*         } */
  /*       } */
  /*       count++; */
  /*     } */
  /*     changed = false; */
  /*   } */

  /*   Console.WriteLine("Number of trial functions = " + count); */

  /*   vector v = new vector(eigens.Count); */
  /*   for (int i = 0; i < eigens.Count; i++) { */
  /*     v[i] = eigens[i]; */
  /*   } */
  /*   return v; */
  /* } */

  // Generate a basis which can be cholesky factorised
  /* public List<matrix> generateBasis(int size) { */
  /*   List<matrix> basis = generateTestFunctions(1); */
  /*   while(basis.Count < size) { */
  /*     List<matrix> tmpbasis = new List<matrix>(basis); */
  /*     tmpbasis.Add(generateTestFunctions(1)[0]); */
  /*     matrix B = generateB(tmpbasis); */
  /*     if (validateB(B)) { */
  /*       basis = tmpbasis; */
  /*     } */
  /*   } */
  /*   return basis; */
  /* } */

  // Ensure the eigenvalues of B are greater than 0.01
  // to prevent too much correlation of basis functions
  /* public bool validateB(matrix B) { */
  /*   vector v = new vector(B.cols); */
  /*   jacobi.eigen(B.copy(),v); */
  /*   if (v[0] > 0.1) { */
  /*     try { */
  /*       new CholeskyDecomposition(B); */
  /*       return true; */
  /*     } catch (CholeskyException e) { */
  /*       Console.WriteLine(e); */
  /*       return false; */
  /*     } */
  /*   } */
  /*   return false; */
  /* } */

  // Generate matrix containing non-linear parameters for the gaussian test
  // functions transformed to center of mass system
  public matrix generateA() {
    matrix A = new matrix(nParticles - 1, nParticles -1);
    List<double> alphas = makeAlphas();
    // Generate values for every entry in the matrix A
    for (int k = 0; k < nParticles - 1; k++) {
      /* for (int l = 0; l < nParticles - 1; l++) { */
      /*   A[k,l] = A_kl(k,l,alphas); */
      /* } */
      for (int l = 0; l < k; l++) {
        double akl = A_kl(k,l,alphas);
        A[l,k] = akl;
        A[k,l] = akl;
      }
      A[k,k] = A_kl(k,k,alphas);
    }

    return A;
  }

  public List<double> makeAlphas() {

    List<double> alphas = new List<double>(); // non-linear parameters
    double min = problem.minGuess(); // Minimum value to look for parameters
    double max = problem.maxGuess(); // Maximum value to look for parameters

    // Prepare randdom values
    for (int i = 0; i < nParticles * (nParticles - 1) / 2; i++) {
      alphas.Add( r.NextDouble() * (max - min) + min);
    }

    return alphas;
  }

  // Generate the value in A belonging to spot A[k,l]
  public double A_kl(int k, int l, List<double> alphas) {
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


  // Generate a list of permutations taking distinguishable particles into
  // account. That is protons and electrons won't be swapped
  public List<Permutation> permutations () {
    // Generate List of List of indistinguishable particles
    List<List<int>> listlist = new List<List<int>>();
    List<int> tmplist = new List<int>();
    Particle p = problem.getParticles()[0];
    tmplist.Add(0);
    for(int i = 1; i < nParticles; i++) {
      if (problem.getParticles()[i].equals(p)) {
        tmplist.Add(i);
      } else {
        listlist.Add(new List<int>(tmplist));
        tmplist = new List<int>();
        tmplist.Add(i);
        p = problem.getParticles()[i];
      }
      if ( i == nParticles - 1) {
        listlist.Add(new List<int>(tmplist));
      }
    }

    // Generate List of List of permutations
    List<List<Permutation>> llp = new List<List<Permutation>>();
    foreach(List<int> list in listlist) {
      List<Permutation> permu = new List<Permutation>();
      permu.Add(new Permutation(list,1,-1)); // Need to change last parameter if boson
      llp.Add(Misc.permutations(permu));
    }

    // List containing all previous permutations
    List<Permutation> lp = new List<Permutation>();

    for (int i = llp.Count - 1; i > -1; i--) {
      List<Permutation> tmp = new List<Permutation>();
      foreach(Permutation per in llp[i]) {
        Permutation clone = per.Clone();
        if (lp.Count == 0) {
          tmp.Add(clone);
        }
        foreach(Permutation peru in lp) {
          Permutation clone2 = per.Clone();
          clone2.Append(peru);
          tmp.Add(clone2);
        }
      }
      lp = new List<Permutation>(tmp);
    }
    return lp;
  }


  // Symmetrize a matrix so that there is anti symmetry on swapping
  // indistinguishable fermions and symmetry on swapping indistinguishable
  // bosons
  public matrix symmetrize(matrix A) {
    // Generate permutations of particles
    List<int> perm = new List<int>();
    for (int i = 0; i < nParticles; i++){
      perm.Add(i);
    }

    // Make a list with 1 permutation. It will be used for getting permutations of
    /* List<Permutation> ll = new List<Permutation>(); */
    /* ll.Add(new Permutation(perm,1,-1)); // fermions */

    matrix result = new matrix(A.rows,A.cols);

    foreach (Permutation ps in permutations()) {
      // Generate C_i according to the permutation and transform to jacobi
      matrix Ci = generateC(ps.ToArray());
      matrix tmp = problem.getUInverse() * Ci * problem.getU();
      /* matrix tmp = problem.getUInverse().transpose() * Ci * problem.getU(); */
      /* matrix tmp = problem.getU().transpose() * Ci * problem.getU(); */
      matrix Pi = new matrix(A.rows, A.cols);

      // Remove last row and column from matrix
      for(int i = 0; i < A.rows; i++) {
        for(int j = 0; j < A.cols; j++) {
         Pi[i,j] = tmp[i,j];
        }
      }

      /* Console.WriteLine("Pi = "); */
      /* Pi.print(); */

      int[] permus = ps.ToArray();
      matrix Tp = new matrix(A.rows,A.cols);
      for(int i = 0; i < Tp.rows; i++) {
        for(int j = 0; j < Tp.cols; j++) {
          double entry = 0;
          matrix U = problem.getU();
          matrix U_inv = problem.getUInverse();
          for(int k = 0; k < permus.Length; k++) {
            entry += U[i,k]*U_inv[permus[k],j];
          }
          Tp[i,j] = entry;
        }
      }
      /* Console.WriteLine("Tp = "); */
      /* Tp.print(); */


      // There might be an issue with size here
      /* result += ps.Parity() * (Pi.transpose() *A * Pi); */
      result += ps.Parity() * (Tp.transpose() * A * Tp.transpose());
    }

    return 1/Math.Sqrt(Misc.factorial(A.rows)) * result;
  }

  // Generate a matrix representing the given permutation
  public matrix generateC(int[] permutation) {
    int length = permutation.Length;
    matrix C = new matrix(length,length);
    for(int i = 0; i < length; i++) {
      for(int j = 0; j < length; j++) {
        if(j == permutation[i]) {
          C[i,j] = 1;
        }
      }
    }
    return C;
  }

  // Generate a number of test functions to work as a basis
  // for the system specified in problem
  public List<matrix> generateTestFunctions(int number) {
    List<matrix> testFunctions = new List<matrix>();
    for (int i=0; i < number; i++) {
      // TODO symmetriseringen virker ikke...
      /* matrix a = symmetrize(generateA()); */
      matrix a = generateA();
      testFunctions.Add(a);
    }

    return testFunctions;
  }

  // Generate the matrix B based on the given test functions
  public matrix generateB(List<matrix> testFunctions) {
    int size = testFunctions.Count;
    matrix B = new matrix(size,size);
    for(int i = 0; i < size; i++) {
      /* for (int j = 0; j < size; j++) { */
      /*   matrix A = testFunctions[i]; */
      /*   matrix C = testFunctions[j]; */
      /*   B[i,j] = mef.overlapElement(A,C); */
      /* } */
      matrix A = testFunctions[i];
      for (int j = 0; j < i; j++) {
        matrix C = testFunctions[j];
        double element = mef.overlapElement(A,C);
        B[i,j] = element;
        B[j,i] = element;
      }
      B[i,i] = mef.overlapElement(A,A);
    }
    return B;
  }

  // Generate the matrix H contain matrix elements of kinetic
  // and potential energy combined
  public matrix generateH(List<matrix> testFunctions) {
    int size = testFunctions.Count;
    matrix H = new matrix(size, size);

    /* for (int i=0; i < size; i++) { */
    /*   for (int j=0; j < size; j++) { */
    /*     H[i,j] = mef.matrixElement(testFunctions[i],testFunctions[j]); */
    /*   } */
    /* } */

    for (int i=0; i < size; i++) {
      for (int j=0; j < i; j++) {
        // H should be symmetric
        double element = mef.matrixElement(testFunctions[i],testFunctions[j]);
        H[i,j] = element;
        H[j,i] = element;
      }
      H[i,i] = mef.matrixElement(testFunctions[i],testFunctions[i]);
    }
    return H;
  }

  // Calculate eigen values of system as specified in problem
  public vector eigenValues(List<matrix> testFunctions) {
    vector v = new vector(testFunctions.Count);
    matrix H = generateH(testFunctions);
    matrix B = generateB(testFunctions);
    matrix L = new CholeskyDecomposition(B).L;
    matrix inv = new QRdecomposition(L).inverse();
    matrix inv_T = new QRdecomposition(L.transpose()).inverse();
    /* matrix i = inv * H; */
    matrix p = inv * H * inv_T;
    jacobi.eigen(p,v);
    return v;
  }
}
