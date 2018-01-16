package pl.intelliseq.genetraps.api.dx.controllers;

import java.io.IOException;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;

import pl.intelliseq.genetraps.api.dx.FilesManager;
import pl.intelliseq.genetraps.api.dx.parser.DxRunner;

@RestController
public class TouchController {

	Logger log = Logger.getLogger(TouchController.class);

	@Autowired
	private FilesManager filesManager;

    @RequestMapping(value = "/touch", method = RequestMethod.GET)
    public String vep() {
    	
    	//return "{\"response\":\"api-vep server is alive\"}";
    	
    	//log.info("Variant: " + variant);
    	
    	try {
			return touch();
		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
			return e.toString();
		}
    }


    @RequestMapping(value = "/mkdir", method = RequestMethod.GET)
	@ResponseBody public String mkDir(
//			@RequestParam(value="owner", required = true) String owner
	){
		return String.format("{\"response\":%s}", filesManager.mkdir());
	}
    
    
    private String touch() throws IOException, InterruptedException {
    	String result = DxRunner.runCommand("dx run touch -iname=test");
    	log.info(result);
    	return this.getJobId(result);
    	
    }

	private String getJobId(String result) {
		try {
			int indexOfJobId = result.indexOf("Job ID") + 8;
			int endIndexOfJobId = result.substring(indexOfJobId).indexOf("\n");
			if (endIndexOfJobId == -1) {
				return result.substring(indexOfJobId);
			}
			return result.substring(indexOfJobId, indexOfJobId + endIndexOfJobId);
		} catch (Exception e) {
			return null;
		}
	}
    
}
