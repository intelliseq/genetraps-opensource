package pl.intelliseq.genetraps.api.security.controllers;

import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;

@RestController
public class CORSController {

    @CrossOrigin(maxAge = 3600)
    @RequestMapping(value = "/", method = RequestMethod.GET)
    public void preflight() {}

}
