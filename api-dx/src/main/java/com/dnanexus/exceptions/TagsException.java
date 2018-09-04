package com.dnanexus.exceptions;

/**
 * A problem with tags (specified in message) has occured
 */
public class TagsException extends RuntimeException {

    public TagsException(String message) {
        super(message);
    }

    @Override
    public String toString() { return "Problem with tags occured: " + this.getMessage(); }

    private static final long serialVersionUID = 2376890343845817594L;
}
