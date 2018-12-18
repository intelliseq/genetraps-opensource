package pl.intelliseq.genetraps.api.dx.controllers;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.web.bind.annotation.*;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;

@RestController
public class PipelineController {

    @Autowired
    private DxApiProcessManager processManager;

    @RequestMapping(value = "/fastqc", method = RequestMethod.POST)
    @ResponseStatus(HttpStatus.ACCEPTED)
    public String fastqc(@RequestParam String fileId) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runFastqc(fileId).getId()).toString();
    }

//    test me not in use
//    @RequestMapping(value = "/bwa", method = RequestMethod.POST, params = {"fastq_file_1", "fastq_file_2"})
//    public String bwa(@RequestParam String fastq_file_1, @RequestParam String fastq_file_2) {
//        return new ObjectMapper().createObjectNode().put("id", processManager.runBwa(fastq_file_1, fastq_file_2).getId()).toString();
//    }

    @RequestMapping(value = "/bwa", method = RequestMethod.POST, params = {"sampleId"})
    @ResponseStatus(HttpStatus.ACCEPTED)
    public String bwa(
            @RequestParam Integer sampleId) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runBwa(sampleId).getId()).toString();
    }

    @RequestMapping(value = "/gatkhc", method = RequestMethod.POST)
    @ResponseStatus(HttpStatus.ACCEPTED)
    public String gatkhc(
            @RequestParam Integer sampleId,
            @RequestParam(required = false) String interval) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runGatkHC(sampleId, interval).getId()).toString();
    }

}
