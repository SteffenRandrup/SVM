public static class Misc {
  public static int factorial(int i) {
    if (i < 0) return 1; // Bad but fuck it
    int result = 1;
    for (int j = 1; j < i+1; j++) {
      result *= j;
    }
    return result;
  }
}
