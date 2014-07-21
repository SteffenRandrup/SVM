using System;
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
    /* particles.Add(new Particle(1,1,3/2)); */
    // Electron
    particles.Add(new Particle(5.485*Math.Pow(10,-4), -1, 1/2));
    particles.Add(new Particle(5.485*Math.Pow(10,-4), -1, 1/2));
    min = 0.1;
    max = 1.5;
    U = generateU();
    U_inv = new QRdecomposition(U).inverse();
    Lambda = generateLambda();

    particles.Sort();
  }

  public ProblemSetup (String filename) {
    // Read values from a file and set them to the ones above
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
