package pl.intelliseq.genetraps.api.dx.controllers;

import org.springframework.security.access.prepost.PreAuthorize;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;

import java.util.UUID;

@RestController
public class WebController {

    @RequestMapping(value = "/foo", method = RequestMethod.GET)
    public String tokenGet() {
        return "read foo " + UUID.randomUUID().toString() +"\n";
    }

    @PreAuthorize("hasAuthority('FOO_WRITE')")
    @RequestMapping(value = "/foo", method = RequestMethod.POST)
    public String tokenPost() {
        return "write foo " + UUID.randomUUID().toString() +"\n";
    }
}
