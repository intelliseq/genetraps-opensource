package pl.intelliseq.genetraps.api.dx.controllers;

import org.springframework.core.io.FileSystemResource;
import org.springframework.core.io.Resource;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.http.ResponseEntity;

@RestController
public class HelloController {
    @RequestMapping(value = "/hello", method = RequestMethod.GET)
    public String hello() {
        return "{\"status\": \"up\"}";
    }

    /* preflight */
    @RequestMapping(
            value = "/**",
            method = RequestMethod.OPTIONS
    )
    public ResponseEntity handle() {
        return new ResponseEntity(HttpStatus.OK);
    }

    @RequestMapping(value = "/api-dx-online-svg-badge", method = RequestMethod.GET)
    public Resource getBadge() {
        Resource resource = new FileSystemResource("src/main/resources/api--dx-online-brightgreen.svg");
        return resource;
    }
}
