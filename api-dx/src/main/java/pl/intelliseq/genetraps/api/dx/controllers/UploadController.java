package pl.intelliseq.genetraps.api.dx.controllers;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;
import pl.intelliseq.genetraps.api.dx.helpers.ProcessManager;
import pl.intelliseq.genetraps.api.dx.models.IseqJSON;

@RestController
public class UploadController {

	private Logger log = Logger.getLogger(UploadController.class);

	@Autowired
    private ProcessManager processManager;

    @RequestMapping(value = "/upload", method = RequestMethod.POST)
    public String upload(
            @RequestParam String url,
            @RequestParam String sampleNumber,
            @RequestParam(required = false) String tag){
        return new IseqJSON("id", processManager.runUrlFetch(url, sampleNumber, tag)).toString();
    }

    @RequestMapping(value = "/describe/{id}", method = RequestMethod.GET)
    public String describe(@PathVariable String id){
        return processManager.runJSONDescribe(id).toString();

    }
    
}
