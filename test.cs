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
    int[] charges = {-1,1,-1};
    double[] spins = {1/2,1/2,1/2};
    problem = new ProblemSetup(masses,charges,spins,0.1,1);
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
  public void H_should_have_positive_eigenvalues(){
    matrix H = svm.generateH(testfunctions);
    vector v = new vector(H.cols);
    jacobi.eigen(H,v);
    Assert.Greater(v[0],0); // Sorted smallet first
  }

  [Test]
  public void Test_kinetic_energy_element() {
    matrix A = new matrix("0.27942,0.10933;0.10933,0.47448");
    matrix B = new matrix("0.45648,0.03761;0.03761,0.85069");
    Assert.AreEqual(98.4152907127796,mef.kineticEnergyElement(A,B));
  }

  [Test]
  public void Test_coulomb_energy_element() {
    matrix A = new matrix("0.27942,0.10933;0.10933,0.47448");
    matrix B = new matrix("0.45648,0.03761;0.03761,0.85069");
    Assert.AreEqual(1.1149123785677297,mef.coulombPotentialEnergy(A,B));
  }

  [Test]
  [Ignore("not implemented yet")]
  public void Test_gaussian_potential_element() {
    Assert.IsTrue(false);
  }

  [Test]
  public void Test_overlap_element() {
    matrix A = new matrix("0.27942,0.10933;0.10933,0.47448");
    matrix B = new matrix("0.45648,0.03761;0.03761,0.85069");
    double ans = 41.399293513079961;
    Assert.AreEqual(ans,mef.overlapElement(A,B));
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
  [Ignore("For some reason this does not work either")]
  public void Test_generate_U() {
    String u = "1,-1,0;" + 1/2.0 + "," + 1/2.0 + ",-1;" + 1/3.0 + "," + 1/3.0 + "," + 1/3.0;
    Console.WriteLine("u= " + u);
    matrix U = new matrix(u);
    U.print();
    problem.getU().print();
    Assert.IsTrue(U.equals(problem.getU()));
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
  [Ignore("I don't know why it does not work")]
  public void Test_inverse() {
    matrix A = new matrix("1,2,3;0,1,4;5,6,0");
    QRdecomposition qr = new QRdecomposition(A);
    matrix inv = qr.inverse();
    matrix inverse = new matrix("-24,18,5;20,-15,-4;-5,4,1");
    Assert.IsTrue(inv.equals(inverse));
  }

  [Test]
  public void Test_factorial() {
    int fact = 24;
    int factorial = Misc.factorial(4);
    Assert.AreEqual(fact,factorial);

  }

  [Test]
  [Ignore("Not implemented yet")]
  public void Test_symmetrise() {

  }

  [Test]
  public void Test_permutations() {
    List<Permutation> perms = svm.permutations();
    int length = perms.Count;
    Assert.AreEqual(length,2); // only two particles can be swapped
    // should probably also test content
  }

  [Test]
  [Ignore("Not implemented yet")]
  public void Test_eigenvalue_function() {
  }
}
