using System;

class MainClass {

  public static void Main (string[] args) {
    StochasticVariationalMethod(args);
  }

  // Program running for stochastic variational method
  private static void StochasticVariationalMethod(string[] args) {
    SVM svm; // Class handling all aspects of SVM
    String filename = "result.txt";
    if (args.Length == 0) {
      /* svm = new SVM(new ProblemSetup("hydrogen.txt")); */
      svm = new SVM(new ProblemSetup());
    } else {
      svm = new SVM(new ProblemSetup(args[0]));
    }
    vector v = svm.run(20);
    v.print();
    String forprint = "";
    for(int i = 0; i < v.size; i++) {
      forprint += v[i] + "\n";
    }

    System.IO.StreamWriter file = new System.IO.StreamWriter(filename);
    file.WriteLine(forprint);
    file.Close();
  }

  // Program running for Newton-Raphson root finding method
  private static void NewtonRaphson () {

  }

}
