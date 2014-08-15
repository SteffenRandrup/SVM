using System;
public class CholeskyException: Exception {

    public CholeskyException() : base() {

    }

    public CholeskyException(String message) : base(message) {

    }

    public CholeskyException(String message, Exception innerException)
      : base(message, innerException) {

    }
}
