package pl.intelliseq.genetraps.api.security.controllers;

import org.springframework.core.io.FileSystemResource;
import org.springframework.core.io.Resource;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;
import javax.servlet.http.HttpServletResponse;
import java.util.Date;

@RestController
public class HelloController {
    @RequestMapping(value = "/hello", method = RequestMethod.GET)
    public String hello() {
        return "{\"status\": \"up\"}";
    }

    @RequestMapping(value = "/api-security-online-svg-badge", method = RequestMethod.GET)
    public Resource getBadge(HttpServletResponse response) {

        Resource resource = new FileSystemResource("src/main/resources/api--security-online-brightgreen.svg");
	response.setHeader("Cache-Control", "no-cache");
	response.setContentType("image/svg+xml");
	Date date= new Date();
	response.setHeader("ETag","" + date.getTime());
        return resource;
    }
}
