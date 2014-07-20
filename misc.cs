using System.Collections.Generic;

public static class Misc {
  public static int factorial(int i) {
    if (i < 0) return 1; // Bad but fuck it
    int result = 1;
    for (int j = 1; j < i+1; j++) {
      result *= j;
    }
    return result;
  }

  // Recursively generate a list of permutations.
  // The permutation class contains the order and parity of a
  // permutation.
  public static List<Permutation> permutations(List<Permutation> permlist) {
    List<Permutation> result = new List<Permutation>();
    foreach (Permutation p in permlist) {
      // If length is 1 there is only one permutation
      if ( p.Length() == 1) {
        result.Add(p.Clone());}
      // If length is 2 swap the content
      else if ( p.Length() == 2) {
        result.Add(p.Clone());
        result.Add(p.swap(0,1));}
      // For every number in permutation remove it and generate
      // permutation of the remaining. Then append them to the number
      else {
        for (int i = 0; i < p.Length(); i++) {
          Permutation pc = p.Clone();
          Permutation pp = pc.Pop(i);
          List<Permutation> list = new List<Permutation>();
          list.Add(pc);
          foreach(Permutation put in permutations(list)) {
            result.Add(pp.Clone().Append(put));
          }
        }
      }
    }
    return result;
  }

}
