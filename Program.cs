using System;

class MainClass
{
  /* private SVM svm; */

  public static void Main (string[] args) {
    /* StochasticVariationalMethod(); */
    SVM svm = new SVM(new ProblemSetup("hydrogen.txt"));
    foreach (Permutation p in svm.perm2()) {
      foreach(int i in p.ToArray()) {
        Console.Write(i);
      }
      Console.WriteLine();
    }
    
  }

  // Program running for stochastic variational method
  /* private static void StochasticVariationalMethod() { */
  /*   if (args.Length == 0) { */
  /*     svm = new SVM(new ProblemSetup()); */
  /*   } else { */
  /*     svm = new SVM(new ProblemSetup(args[0])); */
  /*   } */

  /*   // Get converging eigen values */
  /*   vector eigs = svm.run(10); */
  /*   String forprint = ""; */

  /*   for (int j = 0; j < eigs.size; j++) { */
  /*     forprint += eigs[j] +"\n"; */
  /*   } */

  /*   System.IO.StreamWriter file = new System.IO.StreamWriter(args[0] + ".res"); */
  /*   file.WriteLine(forprint); */
  /*   file.Close(); */
  /* } */

  /* private static void NewtonRaphson () { */
  
  /* } */


}
