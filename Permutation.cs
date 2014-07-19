using System;
using System.Collections.Generic;

public class Permutation{
  
  private int parity = 1;
  private int pf;
  private List<int> perm;

  public Permutation(List<int> perm, int parity, int factor) {
    this.perm = perm;
    this.parity = parity;
    this.pf = factor;
  }

  public int[] ToArray() {
    return perm.ToArray();
  }

  public int Parity() {
    return parity;
  }

  public int Length() {
    return perm.Count;
  }

  public Permutation Clone() {
    return new Permutation(new List<int>(perm), parity, pf);
  }

  public Permutation swap (int i, int j) {
    if ( i != j ) {
      int tmp = perm[i];
      perm[i] = perm[j];
      perm[j] = tmp;
      parity *= pf;
    }
    return this;
  }

  public Permutation Pop(int i) {
    int res = perm[i];
    perm.RemoveAt(i);
    if ( i > 0) {parity *= pf;}
    List<int> l = new List<int>();
    l.Add(res);
    Permutation p = new Permutation(l, parity, pf);
    return p;
  }

  public Permutation Append(Permutation p) {
    perm.AddRange(p.ToArray());
    parity *= p.Parity();
    return this;
  }

  public override String ToString() {
    String s = "";
    foreach(int i in perm) {
      s += i;
    }
    return s;
  }
}
