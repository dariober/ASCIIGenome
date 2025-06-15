package exceptions;

public class InvalidGenomicCoordsException extends Exception {

  public InvalidGenomicCoordsException() {}

  public InvalidGenomicCoordsException(String msg) {
    super(msg);
  }

  /** */
  private static final long serialVersionUID = 1L;
}
