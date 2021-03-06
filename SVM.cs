using System;
using System.Collections.Generic;

public class SVM {

  private ProblemSetup problem; // Class containing all info
  private int nParticles;
  private Random r = new Random();
  private MEF mef;
  List<Permutation> perms;

  public SVM(ProblemSetup problem) {
    this.problem = problem;
    nParticles = problem.getNumberOfParticles();
    mef = new MEF(problem);
    /* perms = permutations(); */
  }

  public List<double> run2(int size, int testsize) {
    double E0 = double.PositiveInfinity;
    List<double> results = new List<double>();
    List<matrix> basis = generateTestFunctions(size);
    bool converged = false;
    int count = 0;
    double difference = double.PositiveInfinity;
    double tolerance = Math.Pow(10,-3);
    /* for (int k = 0; k < 1000; k++) { */
    while(!converged) {
      for(int i = 0; i < size; i++) {
        for (int j = 0; j < testsize; j++) {
          List<matrix> tmpbasis = new List<matrix>(basis);
          tmpbasis[i] = generateA();
          if (validateB(generateB(tmpbasis))) {
            double low = eigenValues(tmpbasis)[0];
            if (low < E0) {
              difference = Math.Abs((E0 - low)/E0);
              E0 = low;
              Console.WriteLine(low);
              basis = tmpbasis;
              results.Add(low);
              if(difference > tolerance) {
                count = 0;
              }
            }
          } else {
            Console.WriteLine("Ignored a testfunction");
          }
          /* results.Add(E0); */
        }
      }
      if (count > 500) {converged = true;}
      count++;
    }
    return results;
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
        if (validateB(generateB(tmpbasis))){
          vector v = eigenValues(tmpbasis);
          if (v[0] < E0) {
            bestBasisSoFar = tmpbasis;
            E0 = v[0];
            Console.WriteLine(E0);
          }
        }
      }
      result[i] = E0;
      basis = bestBasisSoFar;
    }
    return result;
  }

  public matrix permute(matrix M, Permutation p) {
    matrix Ci = generateC(p.ToArray());
    /* matrix tmp = problem.getU().transpose() * Ci * problem.getU(); */
    matrix tmp = problem.getU() * Ci * problem.getUInverse();
    matrix T = new matrix(M.rows, M.cols);
    for(int i = 0; i < M.rows; i++) {
      for(int j = 0; j < M.cols; j++) {
        T[i,j] = tmp[i,j];
      }
    }
    matrix result = T * M * T.transpose();
    return result;
  }


  // Ensure the eigenvalues of B are greater than 0.01
  // to prevent too much correlation of basis functions
  public bool validateB(matrix B) {
    vector v = new vector(B.cols);
    jacobi.eigen(B.copy(),v);
    /* if (v[0] > 0.01) { */
      try {
        new CholeskyDecomposition(B);
        return true;
      } catch (CholeskyException) {
        return false;
      }
    /* } else { */
    /*   B.print(); */
    /*   Console.WriteLine("Eigenvalue " + v[0] + " is not large enough"); */
    /* } */
    /* return false; */
  }

  // Generate matrix containing non-linear parameters for the gaussian test
  // functions transformed to center of mass system
  /* public matrix generateA() { */
  /*   matrix A = new matrix(nParticles - 1, nParticles -1); */
  /*   List<double> alphas = makeAlphas(); */
  /*   // Generate values for every entry in the matrix A */
  /*   for (int k = 0; k < nParticles - 1; k++) { */
  /*     for (int l = 0; l < k; l++) { */
  /*       double akl = A_kl(k,l,alphas); */
  /*       A[l,k] = akl; */
  /*       A[k,l] = akl; */
  /*     } */
  /*     A[k,k] = A_kl(k,k,alphas); */
  /*   } */

  /*   return A; */
  /* } */

  public matrix generateA() {
    matrix A = new matrix(nParticles - 1, nParticles - 1);
    matrix B = new matrix(A.rows,A.cols);
    for (int i = 0; i < A.cols; i++) {
      for (int j = 0; j < A.cols; j++) {
        if (i == j) {
          A[i,j] = Math.Log(1-r.NextDouble())/-1;
        }
        B[i,j] = Math.Log(1-r.NextDouble())/-1;
      }
    }
    matrix Q = new QRdecomposition(B).Q;

    return Q * A * Q.transpose();
  }

  public List<double> makeAlphas() {

    List<double> alphas = new List<double>(); // non-linear parameters
    double min = problem.minGuess(); // Minimum value to look for parameters
    double max = problem.maxGuess(); // Maximum value to look for parameters

    // Prepare randdom values
    for (int i = 0; i < nParticles * (nParticles - 1) / 2; i++) {
      double rateParam = -1 * (max - min);
      double rand = Math.Log(1 - r.NextDouble()) / rateParam;
      /* Console.WriteLine(rand); */
      alphas.Add(rand);
      /* alphas.Add( r.NextDouble() * (max - min) + min); */
    }

    return alphas;
  }

  // Generate the value in A belonging to spot A[k,l]
  public double A_kl(int k, int l, List<double> alphas) {
    double result = 0;


    matrix U_inv = problem.getUInverse(); // Get the inverse U matrix
    List<double> products = new List<double>();


    for (int j = 0; j < U_inv.cols; j++) {
      for (int i = 0; i < j; i++) {
        products.Add((U_inv[i,k] - U_inv[j,k]) * (U_inv[i,l] - U_inv[j,l]));
      }
    }

    for (int i = 0; i < alphas.Count; i++) {
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
      testFunctions.Add(generateA());
    }
    return testFunctions;
  }


  public matrix generateB(List<matrix> testFunctions) {
    int size = testFunctions.Count;
    matrix H = new matrix(size, size);

    for (int i=0; i < size; i++) {
      for (int j=0; j < i; j++) {
        // H should be symmetric
        double element = mef.overlapElement(testFunctions[i],testFunctions[j]);
        H[i,j] = element;
        H[j,i] = element;
      }
      H[i,i] = mef.overlapElement(testFunctions[i],testFunctions[i]);
    }
    return H;
  }

  // Generate B symmetrized
  public matrix generateB2(List<matrix> testFunctions) {
    int size = testFunctions.Count;
    /* List<Permutation> perms = permutations(); */
    /* Console.WriteLine("There are " + perms.Count + " permutations"); */
    matrix B = new matrix(size,size);
    for (int i = 0; i < size; i++) {
      for (int j = 0; j <= i; j++) {
        double element = 0;
        List<matrix> permuted_i = new List<matrix>();
        List<matrix> permuted_j = new List<matrix>();
        List<int> parities = new List<int>();

        foreach(Permutation p in perms) {
          permuted_i.Add(permute(testFunctions[i], p));
          permuted_j.Add(permute(testFunctions[j], p));
          parities.Add(p.Parity());
        }

        for (int k = 0; k < permuted_i.Count; k++) {
          for (int l = 0; l < permuted_j.Count; l++) {
            int factor = perms[k].Parity() * perms[l].Parity();
            element += factor * mef.overlapElement(permuted_i[k],permuted_j[l]);
          }
        }

        B[i,j] = B[j,i] =  element / Math.Sqrt(perms.Count);
      }
    }

    return B;
  }

  // Generate the matrix H contain matrix elements of kinetic
  // and potential energy combined
  public matrix generateH(List<matrix> testFunctions) {
    int size = testFunctions.Count;
    matrix H = new matrix(size, size);

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

  // Generate H symmetrized
  public matrix generateH2(List<matrix> testFunctions) {
    int size = testFunctions.Count;
    /* List<Permutation> perms = permutations(); */
    matrix H = new matrix(size,size);
    for (int i = 0; i < size; i++) {
      for (int j = 0; j <= i; j++) {
        double element = 0;
        List<matrix> permuted_i = new List<matrix>();
        List<matrix> permuted_j = new List<matrix>();
        List<int> parities = new List<int>();

        foreach(Permutation p in perms) {
          permuted_i.Add(permute(testFunctions[i], p));
          permuted_j.Add(permute(testFunctions[j], p));
          parities.Add(p.Parity());
        }

        for (int k = 0; k < permuted_i.Count; k++) {
          for (int l = 0; l < permuted_j.Count; l++) {
            int factor = perms[k].Parity() * perms[l].Parity();
            element += factor * mef.matrixElement(permuted_i[k],permuted_j[l]);
          }
        }

        H[i,j] = H[j,i] = element / Math.Sqrt(perms.Count);
      }
    }

    return H;
  }


  // Calculate eigen values of system as specified in problem
  public vector eigenValues(List<matrix> testFunctions) {
    vector v = new vector(testFunctions.Count);
    /* try { */
      matrix H = generateH(testFunctions);
      matrix B = generateB(testFunctions);
      matrix L = new CholeskyDecomposition(B).L;
      matrix inv = new QRdecomposition(L).inverse();
      matrix inv_T = new QRdecomposition(L.transpose()).inverse();
      matrix p = inv * H * inv_T;
      jacobi.eigen(p,v);
    /* } catch (CholeskyException ce) { */
    /*   Console.WriteLine(ce); */
    /* } */
    return v;
  }

  public void runSymHelium() {
    double E0 = double.PositiveInfinity;
    int size = 30;
    List<matrix> basis = generateTestFunctions(size);
    List<Permutation> p = new List<Permutation>();
    List<int> l1 = new List<int>(); l1.Add(0); l1.Add(1); l1.Add(2);
    List<int> l2 = new List<int>(); l2.Add(1); l2.Add(0); l2.Add(2);
    p.Add(new Permutation(l1, 1, 1));
    p.Add(new Permutation(l2, -1, 1));
    perms = p;

    for (int k = 0; k < 1000; k++) {
      for (int i = 0; i < size; i++) {
        for (int j = 0; j < 100 ; j++) {
          List<matrix> tmpbasis = new List<matrix>(basis);
          tmpbasis[i] = generateA();
          if (validateB(generateB2(tmpbasis))) {
            vector v = new vector(tmpbasis.Count);
            matrix L = new CholeskyDecomposition(generateB2(tmpbasis)).L;
            matrix inv = new QRdecomposition(L).inverse();
            matrix inv_T = new QRdecomposition(L.transpose()).inverse();
            matrix h = generateH2(tmpbasis);
            matrix m = inv * h * inv_T;
            jacobi.eigen(m,v);
            double low = v[0];
            if (low < E0) {
              E0 = low;
              basis = tmpbasis;
              Console.WriteLine(E0);
              if (E0 < -2.2) {
                h.print();
              }
            }
          }
        }
      }
    }

  }
}
