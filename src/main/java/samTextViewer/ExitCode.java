package samTextViewer;

public enum ExitCode {
    CLEAN, // = 0
    CLEAN_NO_FLUSH, // Command exited clean but do not repeat it if user presses <enter>
    NULL, // Set when no command has been executed yet.   
    ERROR,	
}
