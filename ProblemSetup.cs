using System;
using System.Diagnostics;
using System.Collections.Generic;
using math = System.Math;

public class ProblemSetup
{

  private List<Particle> particles = new List<Particle>(); // All the particles go here
  private double min; // Minimum guess value
  private double max; // Maximum guess value
  private matrix U;   // For linear transform r -> x
  private matrix U_inv; // U.inverse
  private matrix Lambda;

  public ProblemSetup () {
    // Proton
    particles.Add(new Particle(1837.47,1,1/2.0));
    // Electron
    particles.Add(new Particle(1.0,-1,1/2.0));
    min = 0;
    max = 20;
    setup();
  }

  public ProblemSetup(double[] masses, int[] charges, double[] spins, double min, double max) {
    Debug.Assert( masses.Length == charges.Length && charges.Length == spins.Length );
    for (int i = 0; i < masses.Length; i++) {
      particles.Add(new Particle(masses[i],charges[i],spins[i]));
    }
    this.min = min;
    this.max = max;
    setup();
  }

  public ProblemSetup (String filename) {
    // Read values from a file and set them to the ones above
    string[] lines = System.IO.File.ReadAllLines(filename);
    List<double> masses = new List<double>();
    List<int> charges = new List<int>();
    List<double> spin = new List<double>();

    foreach(string m in lines[0].Split(',')) {
      masses.Add(double.Parse(m.Trim()));
    }
    foreach(string q in lines[1].Split(',')) {
      charges.Add(int.Parse(q.Trim()));
    }
    if (lines.Length > 2) {
      foreach(string s in lines[2].Split(',')) {
        spin.Add(double.Parse(s.Trim()));
      }
    } else {
      foreach(int i in charges) {
        spin.Add(0);
      }
    }
    if (lines.Length == 4) {
      min = double.Parse(lines[3].Split(',')[0].Trim());
      max = double.Parse(lines[3].Split(',')[1].Trim());
    } else {
      min = 0;
      max = 20;
    }
    if (masses.Count == charges.Count && masses.Count == spin.Count) {
      for(int i = 0; i < masses.Count; i++) {
        particles.Add(new Particle(masses[i],charges[i],spin[i]));
      }
    } else {
      Console.WriteLine("Input is not of equal length -- erros might occour");
    }
    setup();
  }

  private void setup() {
    particles.Sort();
    U = generateU();
    U_inv = new QRdecomposition(U).inverse();
    Lambda = generateLambda();
  }

  // Return all the particles
  public List<Particle> getParticles() {
    return particles;
  }

  // Return the number of particles
  public int getNumberOfParticles () {
    return particles.Count;
  }

  // Minimum guess value for non-linear parameters
  public double minGuess() {
    return min;
  }

  // Maximum guess value for non-linear parameters
  public double maxGuess() {
    return max;
  }

  public matrix getU() {
    return U.copy();
  }

  public matrix getUInverse() {
    return U_inv.copy();
  }

  public matrix getLambda() {
    return Lambda;
  }

  // Generates the U matrix
  private matrix generateU() {
    int L = particles.Count;
    matrix U = new matrix(L, L);
    double massSum = 0;

    for (int i = 0; i < L; i++) {
      massSum += particles[i].getMass();

      for (int j = 0; j <= i; j++) {
        U [i,j] = particles[j].getMass() / massSum;
      }

      if (i + 1 < L) {
        U [i, i + 1] = -1;
      }
    }
    return U;
  }

  // Generates the Lambda matrix
  // TODO: Check lambda
  private matrix generateLambda() {
    int nParticles = particles.Count;
    matrix Lambda = new matrix(nParticles - 1, nParticles - 1);

    for (int i = 0; i < nParticles - 1; i++) {
      for (int j = 0; j < nParticles - 1; j++) {
        double lambda_ij = 0;
        for(int k = 0; k < nParticles; k++) {
          lambda_ij += U[i,k] * U[j,k] / particles[k].getMass();
        }
        Lambda[i,j] = lambda_ij;
      }
    }


    return Lambda;
  }

}
