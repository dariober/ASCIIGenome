package exceptions;

public class InvalidTrackTypeException extends Exception {
    public InvalidTrackTypeException() {

    }
    public InvalidTrackTypeException(String message)
    {
        super(message);
    }
    /** */
    private static final long serialVersionUID = 1L;
}
