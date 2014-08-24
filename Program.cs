using System;
using System.Collections.Generic;

class MainClass {

  public static void Main (string[] args) {
    StochasticVariationalMethod(args);
  }

  // Program running for stochastic variational method
  private static void StochasticVariationalMethod(string[] args) {
    SVM svm; // Class handling all aspects of SVM
    String filename = "result.txt";
    if (args.Length == 0) {
      svm = new SVM(new ProblemSetup());
    } else {
      svm = new SVM(new ProblemSetup(args[0]));
    }
    List<double> v = svm.run2(5,100);
    String forprint = "";
    foreach(double res in v) {
      forprint += res + "\n";
    }

    System.IO.StreamWriter file = new System.IO.StreamWriter(filename);
    file.WriteLine(forprint);
    file.Close();
  }

  // Program running for Newton-Raphson root finding method
  private static void NewtonRaphson () {

  }

}
