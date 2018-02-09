package pl.intelliseq.genetraps.api.dx.controllers;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;
import pl.intelliseq.genetraps.api.dx.helpers.ProcessManager;
import pl.intelliseq.genetraps.api.dx.models.IseqJSON;

import java.util.Arrays;

@RestController
public class FilesController {

	private Logger log = Logger.getLogger(FilesController.class);

	@Autowired
    private ProcessManager processManager;

    @RequestMapping(value = "/upload", method = RequestMethod.POST)
    public String upload(
            @RequestParam String url,
            @RequestParam String sampleNumber,
            @RequestParam(required = false) String... tags){
        log.info(Arrays.toString(tags));
        return new IseqJSON("id", processManager.runUrlFetch(url, sampleNumber, tags)).toString();
    }

    @RequestMapping(value = "/describe/{id}", method = RequestMethod.GET)
    public String describe(@PathVariable String id){
        return processManager.runJSONDescribe(id).toString();

    }

    @RequestMapping(value = "/fastqc", method = RequestMethod.POST)
    public String fastqc(@RequestParam String fileId){
        return new IseqJSON("id", processManager.runFastqc(fileId)).toString();
    }

    @RequestMapping(value = "/bwa", method = RequestMethod.POST)
    public String bwa(
            @RequestParam String left,
            @RequestParam String right,
            @RequestParam String outputFolder
    ){
        return new IseqJSON("id", processManager.runBwa(left, right, outputFolder)).toString();
    }
    
}
