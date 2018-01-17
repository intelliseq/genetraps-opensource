package pl.intelliseq.genetraps.api.dx.controllers;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;
import pl.intelliseq.genetraps.api.dx.FilesManager;
import pl.intelliseq.genetraps.api.dx.models.UrlFetcherJob;
import pl.intelliseq.genetraps.api.dx.models.UrlFetcherState;

@RestController
public class UploadController {

	Logger log = Logger.getLogger(UploadController.class);

	@Autowired
	private FilesManager filesManager;

    @RequestMapping(value = "/upload", method = RequestMethod.POST)
    public String upload(@RequestParam String url, @RequestParam String sampleNumber){
        UrlFetcherJob urlFetcherJob = new UrlFetcherJob(url, sampleNumber);
        return "{\"id\":\""+urlFetcherJob.getId()+"\"}";
    }

    @RequestMapping(value = "/state/{jobId}", method = RequestMethod.GET)
    public String state(@PathVariable String jobId){
        log.info("Controller Job id:"+jobId);
        UrlFetcherState urlFetcherState = new UrlFetcherState(jobId);
        return "{\"state\":\""+urlFetcherState.getState()+"\"}";

    }
    
}
