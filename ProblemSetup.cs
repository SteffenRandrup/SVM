using System;
using System.Collections.Generic;
using math = System.Math;

namespace Projekt
{
	public class ProblemSetup
	{

    private List<Particle> particles = new List<Particle>();

		public ProblemSetup () {
      particles.Add(new Particle(1,1,3/2));
      particles.Add(new Particle(5.485*Math.Pow(10,4), -1, 1/2));
		}

		public ProblemSetup (String filename) {
			// Read values from a file and set them to the ones above
		}
    public List<Particle> getParticles() {
      return particles;
    }

	}
}

