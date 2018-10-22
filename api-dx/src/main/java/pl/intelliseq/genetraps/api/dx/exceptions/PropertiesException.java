package pl.intelliseq.genetraps.api.dx.exceptions;

/**
 * A problem with existence of 'properties' file (specified in message) has occured
 */
public class PropertiesException extends RuntimeException {

    public PropertiesException(String message) {
        super(message);
    }

    @Override
    public String toString() { return "The 'properties' file: " + this.getMessage(); }

    private static final long serialVersionUID = 2376890343328017594L;
}
