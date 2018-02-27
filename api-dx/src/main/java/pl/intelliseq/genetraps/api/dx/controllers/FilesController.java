package pl.intelliseq.genetraps.api.dx.controllers;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;
import pl.intelliseq.genetraps.api.dx.exceptions.DxRunnerException;
import pl.intelliseq.genetraps.api.dx.helpers.ProcessManager;
import pl.intelliseq.genetraps.api.dx.models.IseqJSON;

import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@RestController
public class FilesController {

	private Logger log = Logger.getLogger(FilesController.class);

	private String leftRegex = "(.*?)_1\\.(f(ast)?q(\\.gz)?)";
	private String rightRegex = "(.*?)_2\\.(f(ast)?q(\\.gz)?)";


	@Autowired
    private ProcessManager processManager;


	private String matchFileAndGetName(String filename, String regex) {
        Pattern p = Pattern.compile(regex);
        Matcher m = p.matcher(filename);
        if(m.matches()){
            return m.group(1);
        }
        return null;
    }


	@RequestMapping(value = "/upload-both", method = RequestMethod.POST)
    public String uploadBoth(
        @RequestParam String left,
        @RequestParam String right,
        @RequestParam String sampleNumber,
        @RequestParam(required = false) String... tags){
        String leftName = matchFileAndGetName(left, leftRegex);
        String rightName = matchFileAndGetName(right, rightRegex);
        if(leftName != null && leftName.equals(rightName)){
            String leftId = processManager.runUrlFetch(left, sampleNumber, tags);
            String rightId = processManager.runUrlFetch(right, sampleNumber, tags);

            return new IseqJSON("id", String.format("[\"%s\",\"%s\"]", leftId, rightId)).toString();
        } else {
            throw new DxRunnerException("Incompatible files");
        }
    }

    @RequestMapping(value = "/upload", method = RequestMethod.POST)
    public String upload(
            @RequestParam String url,
            @RequestParam String sampleNumber,
            @RequestParam(required = false) Boolean force,
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
