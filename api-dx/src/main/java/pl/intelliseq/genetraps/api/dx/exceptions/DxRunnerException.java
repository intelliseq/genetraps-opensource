package pl.intelliseq.genetraps.api.dx.exceptions;

@SuppressWarnings("serial")
public class DxRunnerException extends RuntimeException {

    public DxRunnerException(Exception e) {
        super(e);
    }

    public DxRunnerException(String message) {
        super(message);
    }
}
