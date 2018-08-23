package pl.intelliseq.genetraps.api.dx.controllers;

import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
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
            @RequestParam String url,
            @RequestParam String sampleNumber,
            @RequestParam String... tag) {
        log.info(Arrays.toString(tag));

        return new ObjectMapper().createObjectNode().put("id", processManager.runUrlFetch(url, sampleNumber, tag).getId()).toString();
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

    @RequestMapping(value = "/bwa", method = RequestMethod.POST, params = {"fastq_file_1", "fastq_file_2", "reference"})
    public String bwa(@RequestParam String fastq_file_1, @RequestParam String fastq_file_2, @RequestParam(defaultValue = "genetraps-resources:reference/grch38-no-alt/grch38-no-alt.tar") String reference) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runBwa(fastq_file_1, fastq_file_2, reference).getId()).toString();
    }

    @RequestMapping(value = "/bwa", method = RequestMethod.POST, params = {"samples_number", "reference"})
    public String bwa(@RequestParam int samples_number, @RequestParam(defaultValue = "genetraps-resources:reference/grch38-no-alt/grch38-no-alt.tar") String reference) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runBwa(samples_number, reference).getId()).toString();
    }

    @RequestMapping(value = "/mkdir", method = RequestMethod.GET)
    @ResponseBody
    public String mkDir() {
        return String.format("{\"response\":%s}", filesManager.mkdir());
    }

    @RequestMapping(value = "/sample/{no}/ls", method = RequestMethod.GET)
    public String samplels(@PathVariable("no") int samples_number) {
        return processManager.sampleLs(samples_number);
    }

    @RequestMapping(value = "/sample/{no}/revls", method = RequestMethod.GET)
    public String samplerevls(@PathVariable("no") int samples_number) {
        return processManager.sampleRevLs(samples_number);
    }
}
