package pl.intelliseq.genetraps.api.dx.controllers;

import org.springframework.core.io.FileSystemResource;
import org.springframework.core.io.Resource;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.web.servlet.view.RedirectView;

import javax.servlet.http.HttpServletResponse;
import java.text.SimpleDateFormat;
import java.util.Date;

@RestController
public class HelloController {

    @RequestMapping(value = "/hello", method = RequestMethod.GET)
    public String hello() {
        return "{\"status\": \"up\"}";
    }

    @RequestMapping(value = "/secure-hello", method = RequestMethod.GET)
    public String secureHello() {
        return "{\"status\": \"ok\"}";
    }

    /* preflight */
//    @RequestMapping(
//            value = "/**",
//            method = RequestMethod.OPTIONS
//    )
//    public ResponseEntity handle() {
//        return new ResponseEntity(HttpStatus.OK);
//    }

    @RequestMapping(value = "/", method = RequestMethod.GET)
    public RedirectView root() {
        return new RedirectView("/swagger-ui.html");
    }

    @RequestMapping(value = "/status", method = RequestMethod.GET, produces = "image/svg+xml")
    public Resource getBadge(final HttpServletResponse response) {
        SimpleDateFormat format = new SimpleDateFormat("EEE, dd MMM yyyy HH:mm:ss zzz");
        response.setHeader("Cache-Control", "no-cache");
        response.setHeader("Last-Modified", format.format(new Date()));
        Resource resource = new FileSystemResource("src/main/resources/api--dx-online-brightgreen.svg");
        return resource;
    }

//    @RequestMapping(value = "/swagger", method = RequestMethod.GET)
//    public Resource getSwagger() {
//        Resource resource = new FileSystemResource("src/main/resources/swagger/api-doc.yaml");
//        return resource;
//    }
}
