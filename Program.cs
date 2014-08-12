using System;

class MainClass {

  public static void Main (string[] args) {
    StochasticVariationalMethod(args);
  }

  // Program running for stochastic variational method
  private static void StochasticVariationalMethod(string[] args) {
    SVM svm; // Class handling all aspects of SVM
    if (args.Length == 0) {
      svm = new SVM(new ProblemSetup("hydrogen.txt"));
    } else {
      svm = new SVM(new ProblemSetup(args[0]));
    }
    for (int i = 0; i < 10; i++) {
      vector eigs = svm.run(7);
    }

    String forprint = "";

    for (int j = 0; j < eigs.size; j++) {
      forprint += eigs[j] +"\n";
    }

    /* System.IO.StreamWriter file = new System.IO.StreamWriter(args[0] + ".res"); */
    System.IO.StreamWriter file = new System.IO.StreamWriter("hydrogen.txt.res");
    file.WriteLine(forprint);
    file.Close();
  }

  // Program running for Newton-Raphson root finding method
  private static void NewtonRaphson () {

  }

}
