using System;

public class Particle : IComparable {
  private double mass;
  private int charge;
  private double spin;

  public Particle (double mass, int charge, double spin) {
   this.mass = mass;
   this.charge = charge;
   this.spin = spin;
  }

  public int CompareTo(Object o) {
    Particle p = o as Particle;
    if (p != null) {
      if (p.getMass() == mass) {
        if (p.getCharge() == charge) {
          if (p.getSpin() == spin) {
            return 0;
          }
          return spin.CompareTo(p.getSpin());
        }
        return charge.CompareTo(p.getCharge());
      }
      return mass.CompareTo(p.getMass());
    } else
        throw new ArgumentException("Object is not a Particle");
  }

  public double getMass() {
    return mass;
  }
  public int getCharge() {
    return charge;
  }
  public double getSpin() {
    return spin;
  }

  public bool equals(Particle p) {
    if (CompareTo(p) == 0) {
      return true;
    }
    return false;
  }
}
