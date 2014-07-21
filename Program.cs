using System;

class MainClass {

  public static void Main (string[] args) {
    StochasticVariationalMethod(args);
  }

  // Program running for stochastic variational method
  private static void StochasticVariationalMethod(string[] args) {
    SVM svm; // Class handling all aspects of SVM
    if (args.Length == 0) {
      svm = new SVM(new ProblemSetup());
    } else {
      svm = new SVM(new ProblemSetup(args[0]));
    }

    // Get converging eigen values
    vector eigs = svm.run(5);
    eigs.print();
    /* String forprint = ""; */

    /* for (int j = 0; j < eigs.size; j++) { */
    /*   forprint += eigs[j] +"\n"; */
    /* } */

    /* System.IO.StreamWriter file = new System.IO.StreamWriter(args[0] + ".res"); */
    /* file.WriteLine(forprint); */
    /* file.Close(); */
  }

  // Program running for Newton-Raphson root finding method
  private static void NewtonRaphson () {

  }

}
