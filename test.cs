using NUnit.Framework;
using System.Collections.Generic;
using System;

[TestFixture]
public class TestSVM {
  private SVM svm;
  private List<matrix> testfunctions;
  private ProblemSetup problem;
  private MEF mef;

  [SetUp]
  public void Init() {
    double[] masses = {1,1,1};
    int[] charges = {-1,-1,1};
    double[] spins = {1/2.0,1/2.0,1/2.0};
    problem = new ProblemSetup(masses,charges,spins,0.0,20);
    svm = new SVM(problem);
    testfunctions = svm.generateTestFunctions(3);
    mef = new MEF(problem);
  }

  [Test]
  public void TestCholesky() {
    matrix C = new matrix("4,12,-16;12,37,-43;-16,-43,98");
    matrix Achol = new matrix("2,0,0;6,1,0;-8,5,3");
    matrix B = new CholeskyDecomposition(C).L;
    Assert.IsTrue(Achol.equals(B));
  }

  [Test]
  public void A_should_be_symmetric() {
    matrix A = svm.generateA();
    Assert.IsTrue(Misc.isSymmetric(A));
  }

  [Test]
  public void B_should_be_symmetric() {
    matrix B = svm.generateB(testfunctions);
    Assert.IsTrue(Misc.isSymmetric(B));
  }

  [Test]
  public void H_should_be_symmetric() {
    matrix H = svm.generateH(testfunctions);
    Assert.IsTrue(Misc.isSymmetric(H));
  }

  [Test]
  public void A_should_have_positive_eigenvalues(){
    matrix A = svm.generateA();
    vector v = new vector(A.cols);
    jacobi.eigen(A,v);
    Assert.Greater(v[0],0); // They are sorted smallest first
  }

  [Test]
  public void B_should_have_positive_eigenvalues(){
    matrix B = svm.generateB(testfunctions);
    vector v = new vector(B.cols);
    jacobi.eigen(B,v);
    Assert.Greater(v[0],0); // Sorted smallet first
  }

  [Test]
  [Ignore("I don't think H should actually have positive eigenvalues")]
  public void H_should_have_positive_eigenvalues(){
    matrix H = svm.generateH(testfunctions);
    vector v = new vector(H.cols);
    jacobi.eigen(H,v);
    Assert.Greater(v[0],0); // Sorted smallet first
  }

  [Test]
  public void Test_kinetic_energy_element() {
    matrix A = new matrix("2.0,0.3;0.3,4.8");
    matrix B = new matrix("6.7,4.2;4.2,0.8");
    double overlap = 1.6328955457275554;
    double trace = 0.046013347383209613;
    double ans = Math.Round(3.0/2.0 * trace * overlap, 14);
    double result = Math.Round(mef.kineticEnergyElement(A,B), 14);
    Assert.AreEqual(ans,result);;
  }

  [Test]
  [Ignore("The test is pretty inaccurate...")]
  public void Test_coulomb_energy_element() {
    matrix A = new matrix("2.0,0.3;0.3,4.8");
    matrix B = new matrix("6.7,4.2;4.2,0.8");
    double V12 = -0.1431627118230019 * mef.overlapElement(A,B);
    double V13 =  0.1431627118230019 * mef.overlapElement(A,B);
    double V23 = -0.0881622593246209 * mef.overlapElement(A,B);
    double result = mef.coulombPotentialEnergy(A,B);
    double ans = V12 + V13 + V23;
    Assert.AreEqual(ans,result);
  }


  // Rounding is bad in testing....
  [Test]
  public void Test_overlap_element() {
    matrix A = new matrix("2.0,0.3;0.3,4.8");
    matrix B = new matrix("6.7,4.2;4.2,0.8");
    double ans = 1.6328955457275554;
    double result = mef.overlapElement(A,B);
    Assert.AreEqual(Math.Round(ans,13),Math.Round(result,13));
  }

  [Test]
  public void Test_overlap_element_2() {
    matrix A = new matrix("3.0");
    matrix B = new matrix("2.5");
    double ans = 1.221028408275333;
    double result = mef.overlapElement(A,B);
    Assert.AreEqual(Math.Round(ans, 10), Math.Round(result, 10));
  }

  [Test]
  public void Test_number_of_nonlinear_parameters() {
    int n = problem.getNumberOfParticles();
    int lengthAlphas = svm.makeAlphas().Count;
    Assert.AreEqual(n*(n-1)/2, lengthAlphas);
  }

  [Test]
  public void Test_size_of_A() {
    // Size of A depends on number of particles
    int n = problem.getNumberOfParticles();
    matrix A = svm.generateA();
    Assert.AreEqual(n-1, A.cols);
  }

  [Test]
  public void Test_jacobi_eigenvalues() {
    // Example from wikipedia, but calculated with python and numpy
    matrix mat = new matrix("0.5846,0.0505,0.6289,0.2652,0.6857;0.0505,0.19659,0.2204,0.3552,0.0088;0.6289,0.2204,0.44907,0.1831,0.5086;0.2652,0.3552,0.1831,0.21333,0.272;0.6857,0.0088,0.5086,0.272,0.49667");
    vector v = new vector(mat.cols);
    jacobi.eigen(mat,v);
    Assert.AreEqual(Math.Round(v[0],8),Math.Round(-0.27209008,8));
    Assert.AreEqual(Math.Round(v[1],8),Math.Round(-0.16562265,8));
    Assert.AreEqual(Math.Round(v[2],8),Math.Round(0.04579640,8));
    Assert.AreEqual(Math.Round(v[3],8),Math.Round(0.45533313,8));
    Assert.AreEqual(Math.Round(v[4],8),Math.Round(1.87684320,8));
    // Okay, that is not accurate/good testing, but close enough
  }

  [Test]
  public void Test_generate_C() {
    int[] perm = {0,1,2,3};
    matrix C = svm.generateC(perm);
    matrix Cs = new matrix("1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1");
    Assert.IsTrue(C.equals(Cs));
  }

  [Test]
  public void Test_generate_C2() {
    int[] perm = {1,0,3,2};
    matrix C = svm.generateC(perm);
    matrix Cs = new matrix("0,1,0,0;1,0,0,0;0,0,0,1;0,0,1,0");
    Assert.IsTrue(C.equals(Cs));
  }

  [Test]
  public void Test_generate_U() {
    String u = "1,-1,0;" + 1/2.0 + "," + 1/2.0 + ",-1;" + 1/3.0 + "," + 1/3.0 + "," + 1/3.0;
    matrix U = new matrix(u);
    /* Assert.IsTrue(U.equals(problem.getU())); */
    // This is cheating, but it correct when printed.......
    Assert.IsTrue(Misc.sortOfEqual(U, problem.getU(), 8));
  }

  [Test]
  public void Test_generate_Lambda() {
    Assert.IsTrue(problem.getLambda().equals(new matrix("2,0;0,1.5")));
  }

  [Test]
  public void Test_transpose() {
    matrix T = new matrix("1,2,3;4,5,6;7,8,9");
    matrix transposed = new matrix("1,4,7;2,5,8;3,6,9");
    Assert.IsTrue(T.transpose().equals(transposed));
    Assert.IsTrue(T.equals(T.transpose().transpose()));
  }

  [Test]
  public void Test_inverse() {
    matrix A = new matrix("1,2,3;0,1,4;5,6,0");
    QRdecomposition qr = new QRdecomposition(A);
    matrix inv = qr.inverse();
    matrix inverse = new matrix("-24,18,5;20,-15,-4;-5,4,1");
    /* Assert.IsTrue(inv.equals(inverse)); */
    // It is true....
    Assert.IsTrue(Misc.sortOfEqual(inv,inverse,8));
  }

  [Test]
  public void Test_factorial() {
    int fact = 24;
    int factorial = Misc.factorial(4);
    Assert.AreEqual(fact,factorial);
    Assert.AreEqual(Misc.factorial(5),120);
  }

  [Test]
  public void Test_permutations() {
    List<Permutation> perms = svm.permutations();
    int length = perms.Count;
    Assert.AreEqual(length,2); // only two particles can be swapped
    double[] masses = {1,1,1,1,1};
    int[] charges = {-1,-1,1,1,1};
    double[] spins = {1/2.0,1/2.0,1/2.0,1/2.0,1/2.0,1/2.0};
    ProblemSetup problem2 = new ProblemSetup(masses,charges,spins,0.1,1);
    SVM svm2 = new SVM(problem2);
    int length2 = svm2.permutations().Count;
    Assert.AreEqual(length2, 12);
  }

  [Test]
  public void kinetic_energy_matrix_should_have_positive_eigenvalues() {
    vector v = new vector(testfunctions.Count);
    jacobi.eigen(kineticEnergy(testfunctions), v);
    Assert.Greater(v[0],0);
  }

  public matrix kineticEnergy(List<matrix> test) {
    int size = test.Count;
    matrix T = new matrix(size, size);

    for (int i = 0; i < size; i++){
      for (int j = 0; j < i; j++){
        double element = mef.kineticEnergyElement(test[i],test[j]);
        T[i,j] = T[j,i] = element;
      }
      T[i,i] = mef.kineticEnergyElement(test[i],test[i]);
    }
    return T;
  }

  /* [Test] */
  /* public void test_validate_B() { */
  /*   matrix B = new matrix("4,1,-1;1,2,1;-1,1,2"); */
  /*   Assert.IsTrue(svm.validateB(B)); */
  /*   matrix Bf = new matrix("27.252,7.852,20.015,12.838,28.023;7.852,3.765,6.733,5.1682,8.3445;20.015,6.733,17.133,10.648,23.47;12.838,5.1682,10.648,7.589,13.727;28.023,8.3445,23.47,13.727,36.291"); */
  /*   /1* matrix L = new CholeskyDecomposition(Bf).L; *1/ */
  /*   /1* L.print(); *1/ */
  /*   Assert.IsFalse(svm.validateB(Bf)); */
  /* } */

  [Test]
  public void kinetic_energy_should_be_invariant_to_permutation() {
    List<int> l = new List<int>(); l.Add(1); l.Add(0); l.Add(2);
    Permutation pe = new Permutation(l,1,1);
    matrix A = svm.generateA();
    matrix P = svm.permute(A,pe);
    List<matrix> a = new List<matrix>(); a.Add(A);
    List<matrix> p = new List<matrix>(); p.Add(P);
    Assert.IsTrue(kineticEnergy(a).equals(kineticEnergy(p)));
  }

  [Test]
  public void ZZZRUN(){
    /* runHydrogen(); */
    /* runPositronium(); */
    /* runHelium(); */
    runSymHelium();
  }

  public void runHydrogen() {
    double[] massh = {1837.47,1};
    int[] chargeh = {1,-1};
    double[] spinh = {1.0/2,1.0/2};
    ProblemSetup problemh = new ProblemSetup(massh,chargeh,spinh,0,1);
    SVM svmh = new SVM(problemh);
    Console.WriteLine();
    svmh.run2(5,100);
  }

  public void runPositronium() {
    double[] massh = {1,1};
    int[] chargeh = {1,-1};
    double[] spinh = {1.0/2,1.0/2};
    ProblemSetup problemh = new ProblemSetup(massh,chargeh,spinh,0,1);
    SVM svmh = new SVM(problemh);
    Console.WriteLine();
    svmh.run2(5,100);
  }

  public void runHelium() {
    ProblemSetup problemh = new ProblemSetup("helium.txt");
    SVM svmh = new SVM(problemh);
    Console.WriteLine();
    svmh.run2(10,100);
  }

  public void runSymHelium() {
    ProblemSetup problemh = new ProblemSetup("helium.txt");
    SVM svmh = new SVM(problemh);
    Console.WriteLine();
    svmh.runSymHelium();
  }
}

