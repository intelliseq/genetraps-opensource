package pl.intelliseq.genetraps.api.dx.controllers;

import org.apache.log4j.Logger;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.RestController;

@RestController
public class HelloController {

	private Logger log = Logger.getLogger(HelloController.class);
	
    @RequestMapping(value = "/hello", method = RequestMethod.GET)
    public String hello() {
        return "{\"status\":\"up\"}";
    }

    @RequestMapping(value = "/hello", method = RequestMethod.POST)
	public String helloPost(@RequestParam String name){
    	return "Hello "+name;
	}
    
}
