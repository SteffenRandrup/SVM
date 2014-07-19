using System;

class MainClass
{

  public static void Main (string[] args) {
    SVM svm = new SVM(new ProblemSetup());
    svm.generateA().print();
  }


}
