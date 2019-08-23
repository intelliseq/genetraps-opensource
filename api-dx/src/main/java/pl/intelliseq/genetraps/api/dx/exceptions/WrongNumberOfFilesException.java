package pl.intelliseq.genetraps.api.dx.exceptions;

/**
 * The number of files found is incorrect/unexpected
 */
public class WrongNumberOfFilesException extends RuntimeException {

    public WrongNumberOfFilesException(String message) {
        super(message);
    }

    @Override
    public String toString() {
        return "WrongNumberOfFiles: " + this.getMessage();
    }

    private static final long serialVersionUID = 2676890341627817594L;
}
