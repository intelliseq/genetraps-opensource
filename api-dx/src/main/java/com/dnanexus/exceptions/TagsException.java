package com.dnanexus.exceptions;

/**
 * The number of files found is incorrect/unexpected
 */
public class TagsException extends RuntimeException {

    public TagsException(String message) {
        super(message);
    }

    @Override
    public String toString() { return "Problem with tags occured: " + this.getMessage(); }

    private static final long serialVersionUID = 2376890343845817594L;
}
