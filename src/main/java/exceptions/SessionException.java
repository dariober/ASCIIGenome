package exceptions;

public class SessionException extends Exception {
    public SessionException() {

    }
    public SessionException(String message)
    {
        super(message);
    }
    private static final long serialVersionUID = 1L;
}
