package pl.intelliseq.genetraps.api.dx;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.concurrent.Executors;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;

import pl.intelliseq.genetraps.api.dx.parser.DxRunner;

@RestController
public class TouchController {

	Logger log = Logger.getLogger(TouchController.class);

	@Autowired
	private FilesManager filesManager;
	
    @RequestMapping("/hello")
    public String hello() {
        return "{\"status\":\"up\"}";
    }
    
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
