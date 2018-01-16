package pl.intelliseq.genetraps.api.dx.controllers;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;
import pl.intelliseq.genetraps.api.dx.FilesManager;
import pl.intelliseq.genetraps.api.dx.parser.DxRunner;

import java.io.IOException;

@RestController
public class HelloController {

	Logger log = Logger.getLogger(HelloController.class);

	@Autowired
	private FilesManager filesManager;
	
    @RequestMapping(value = "/hello", method = RequestMethod.GET)
    public String hello() {
        return "{\"status\":\"up\"}";
    }

    @RequestMapping(value = "/hello", method = RequestMethod.POST)
	public String helloPost(@RequestParam String name){
    	return "Hello "+name;
	}
    
}
