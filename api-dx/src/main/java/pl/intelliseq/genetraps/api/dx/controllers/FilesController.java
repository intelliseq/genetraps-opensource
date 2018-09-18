package pl.intelliseq.genetraps.api.dx.controllers;

import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.oauth2.provider.OAuth2Authentication;
import org.springframework.web.bind.annotation.*;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;

import java.util.Arrays;

@RestController
@Log4j2
public class FilesController {

    @Autowired
    private DxApiProcessManager processManager;

    @Autowired
    private FilesManager filesManager;

//TODO: For future features
//    private String matchFileAndGetName(String filename, String regex) {
//        Pattern p = Pattern.compile(regex);
//        Matcher m = p.matcher(filename);
//        if (m.matches()) {
//            return m.group(1);
//        }
//        return null;
//    }
//
//
//    @RequestMapping(value = "/upload-both", method = RequestMethod.POST)
//    public String uploadBoth(
//            @RequestParam String left,
//            @RequestParam String right,
//            @RequestParam String sampleNumber,
//            @RequestParam(required = false) String... tag) {
//
//        String leftRegex = "(.*?)_1\\.(f(ast)?q(\\.gz)?)";
//        String leftName = matchFileAndGetName(left, leftRegex);
//
//        String rightRegex = "(.*?)_2\\.(f(ast)?q(\\.gz)?)";
//        String rightName = matchFileAndGetName(right, rightRegex);
//
//        if (leftName != null && leftName.equals(rightName)) {
//            DXJob leftId = processManager.runUrlFetch(left, sampleNumber, tag);
//            DXJob rightId = processManager.runUrlFetch(right, sampleNumber, tag);
//
//            return String.format("[\"%s\",\"%s\"]", leftId.getId(), rightId.getId());
//        } else {
//            throw new DxRunnerException("Incompatible files");
//        }
//    }

    @RequestMapping(value = "/upload", method = RequestMethod.POST)
    public String upload(
            OAuth2Authentication auth,
            @RequestParam String url,
            @RequestParam String sampleid,
            @RequestParam String... tag) {
        log.info("upload");
        log.debug(Arrays.toString(tag));
//        log.debug(auth.getUserAuthentication().getPrincipal().toString());

        return new ObjectMapper().createObjectNode().put("id", processManager.runUrlFetch(url, sampleid, tag).getId()).toString();
    }

    @RequestMapping(value = "/describe/{id}", method = RequestMethod.GET)
    public String describe(
            @PathVariable String id) {
        return processManager.JSONDescribe(id).toString();
    }

    @RequestMapping(value = "/fastqc", method = RequestMethod.POST)
    public String fastqc(@RequestParam String fileId) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runFastqc(fileId).getId()).toString();
    }

    @RequestMapping(value = "/bwa", method = RequestMethod.POST, params = {"fastq_file_1", "fastq_file_2"})
    public String bwa(@RequestParam String fastq_file_1, @RequestParam String fastq_file_2) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runBwa(fastq_file_1, fastq_file_2).getId()).toString();
    }

    @RequestMapping(value = "/bwa", method = RequestMethod.POST, params = {"sampleid"})
    public String bwa(@RequestParam int sampleid) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runBwa(sampleid).getId()).toString();
    }

    @RequestMapping(value = "/gatkhc", method = RequestMethod.POST)
    public String gatkhc(@RequestParam int sampleid,
                         @RequestParam(required = false) String interval) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runGatkHC(sampleid, interval).getId()).toString();
    }

    @RequestMapping(value = "/mkdir", method = RequestMethod.GET)
    @ResponseBody
    public String mkDir() {
        return String.format("{\"response\":%s}", filesManager.mkdir());
    }

    @RequestMapping(value = "/sample/{no}/ls", method = RequestMethod.GET)
    public String samplels(@PathVariable("no") int sampleid) {
        return processManager.sampleLs(sampleid);
    }

    @RequestMapping(value = "/sample/{no}/revls", method = RequestMethod.GET)
    public String samplerevls(@PathVariable("no") int sampleid) {
        //Reverse ls - keys swapped with values
        return processManager.sampleRevLs(sampleid);
    }
}
