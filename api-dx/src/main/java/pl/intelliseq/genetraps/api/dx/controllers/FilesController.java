package pl.intelliseq.genetraps.api.dx.controllers;

import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.oauth2.provider.OAuth2Authentication;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import pl.intelliseq.genetraps.api.dx.Roles;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

@RestController
@Log4j2
public class FilesController {

    @Autowired
    private DxApiProcessManager processManager;

    @Autowired
    private FilesManager filesManager;

    @Autowired
    private AuroraDBManager auroraDBManager;

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

    @RequestMapping(value = "/mkdir", method = RequestMethod.GET)
    @ResponseBody
    public String mkDir(OAuth2Authentication auth) {
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        Integer sampleId = filesManager.mkdir();

        auroraDBManager.setUserPrivilegesToSample(userId, sampleId, Roles.ADMIN);

        return String.format("{\"response\":%s}", sampleId);
    }

    @RequestMapping(value = "/upload", method = RequestMethod.POST)
    public String upload(
            OAuth2Authentication auth,
            @RequestParam String url,
            @RequestParam String sampleid,
            @RequestParam String... tag) {
        log.info("upload");
        log.debug(Arrays.toString(tag));
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());
        log.debug(userId);

        return new ObjectMapper().createObjectNode().put("id", processManager.runUrlFetch(url, sampleid, tag).getId()).toString();
    }

    @RequestMapping(value = "/uploadfile", method = RequestMethod.POST)
    public String uploadfile(
            @RequestParam("file") MultipartFile file,
            @RequestParam(value = "sampleid") int sampleid,
            @RequestParam(value = "newfilename", required = false) String newfilename,
            @RequestParam(value = "tag", required = false) List<String> tags) {

        try {
            return new ObjectMapper().createObjectNode().put("id", processManager.runUploadFile(file, sampleid, newfilename, tags)).toString();
        } catch (IOException e) {
            return new ObjectMapper().createObjectNode().put("id", "").toString();
        }
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

//    test me not in use
//    @RequestMapping(value = "/bwa", method = RequestMethod.POST, params = {"fastq_file_1", "fastq_file_2"})
//    public String bwa(@RequestParam String fastq_file_1, @RequestParam String fastq_file_2) {
//        return new ObjectMapper().createObjectNode().put("id", processManager.runBwa(fastq_file_1, fastq_file_2).getId()).toString();
//    }

    @RequestMapping(value = "/bwa", method = RequestMethod.POST, params = {"sampleid"})
    public String bwa(@RequestParam int sampleid) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runBwa(sampleid).getId()).toString();
    }

    @RequestMapping(value = "/gatkhc", method = RequestMethod.POST)
    public String gatkhc(@RequestParam int sampleid,
                         @RequestParam(required = false) String interval) {
        return new ObjectMapper().createObjectNode().put("id", processManager.runGatkHC(sampleid, interval).getId()).toString();
    }

    @RequestMapping(value = "/sample/{no}/ls", method = RequestMethod.GET)
    public String samplels(@PathVariable("no") int sampleid,
                            @RequestParam(required = false, defaultValue = "false") boolean byNames) {
        return processManager.sampleLs(sampleid, byNames);
    }

}
