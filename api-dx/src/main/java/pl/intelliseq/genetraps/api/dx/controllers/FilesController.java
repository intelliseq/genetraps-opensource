package pl.intelliseq.genetraps.api.dx.controllers;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.log4j.Log4j2;
import org.json.JSONObject;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.security.oauth2.provider.OAuth2Authentication;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import pl.intelliseq.genetraps.api.dx.exceptions.PropertiesException;
import pl.intelliseq.genetraps.api.dx.helpers.AWSApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;
import springfox.documentation.annotations.ApiIgnore;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

@RestController
@Log4j2
public class FilesController {

    @Autowired
    private DxApiProcessManager dxApiProcessManager;

    @Autowired
    private AWSApiProcessManager awsApiProcessManager;

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
//            DXJob leftId = dxApiProcessManager.runUrlFetch(left, sampleNumber, tag);
//            DXJob rightId = dxApiProcessManager.runUrlFetch(right, sampleNumber, tag);
//
//            return String.format("[\"%s\",\"%s\"]", leftId.getId(), rightId.getId());
//        } else {
//            throw new DxRunnerException("Incompatible files");
//        }
//    }

    // AWS S3
    @RequestMapping(value = "/sample/create", method = RequestMethod.GET)
    @ResponseBody
    public String createFolder(@ApiIgnore OAuth2Authentication auth) {
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        Integer sampleId = filesManager.mkdir();

        //TO DO make creating samples synchronized with DB i.e. setUserPrivilegesToSample and removing it
//        auroraDBManager.setUserPrivilegesToSample(userId, sampleId, Roles.ADMIN);

        return String.format("{\"response\":\"%s\"}", sampleId);
    }

    // AWS S3
    @RequestMapping(value = "/wdl", method = RequestMethod.POST)
    @ResponseStatus(HttpStatus.ACCEPTED)
    @ResponseBody
    public String wdl(
            @ApiIgnore OAuth2Authentication auth,
            @RequestParam String workflowUrl,
            @RequestParam JSONObject workflowInputs,
            @RequestParam JSONObject labels,
            @RequestParam(name = "req-out", required = false, defaultValue = "{}") JSONObject requestedOutputs) {
        try {
//            if(requestedOutputs.length() == 0)
//                throw new Exception("No requested output set");
            Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());
            return new ObjectMapper().createObjectNode().put("id", awsApiProcessManager.runWdl(userId, workflowUrl, workflowInputs, labels, requestedOutputs)).toString();
        } catch (Exception e) {
            return new ObjectMapper().createObjectNode().put("id", e.getMessage()).toString();
        }
    }

    @RequestMapping(value = "/sample/{id}/urlupload", method = RequestMethod.POST)
    @ResponseStatus(HttpStatus.ACCEPTED)
    public String upload(
            @ApiIgnore OAuth2Authentication auth,
            @PathVariable Integer id,
            @RequestParam String url,
            @RequestParam String... tag) {
        log.info("upload");
        log.debug(Arrays.toString(tag));
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());
        log.debug(userId);

        return new ObjectMapper().createObjectNode().put("id", dxApiProcessManager.runUrlFetch(url, id, tag).getId()).toString();
    }

    // AWS S3
    @RequestMapping(value = "/sample/{id}/file/upload", method = RequestMethod.POST)
    @ResponseStatus(HttpStatus.ACCEPTED)
    public String uploadfile(
            @PathVariable Integer id,
            @RequestParam MultipartFile file,
            @RequestParam(value = "new-name", required = false) String fileName,
            @RequestParam(required = false) List<String> tag) {

        try {
            return new ObjectMapper().createObjectNode().put("id", awsApiProcessManager.runFileUpload(file, id, fileName, tag)).toString();
        } catch (Exception e) {
            return new ObjectMapper().createObjectNode().put("id", e.getMessage()).toString();
        }
    }

    // AWS S3
    @RequestMapping(value = "/sample/{id}/ls", method = RequestMethod.GET)
    public String sampleLs(
            @PathVariable Integer id,
//            @RequestParam(required = false, defaultValue = "false") boolean byNames,
            @RequestParam(required = false, defaultValue = "")  String dir) {
        try {
            return awsApiProcessManager.runSampleLs(id, dir).toString();
        } catch (InterruptedException e) {
            // if error, returns err message with key: /error
            // normal keys doesn't start with '/' at the beginning
            return String.format("{\"/error\":\"%s\"}", e.getMessage());
        }
    }

    //TODO: is there an equivalent in aws??
//    @RequestMapping(value = "/sample/{id}/describe", method = RequestMethod.GET)
//    public String describe(
//            @PathVariable Integer id) {
//        return dxApiProcessManager.JSONDescribe(id).toString();
//    }

    @RequestMapping(value = "/sample/{id}/properties", method = RequestMethod.POST)
    @ResponseStatus(HttpStatus.CREATED)
    public String samplePropertiesPost(
            @PathVariable Integer id,
            @RequestBody LinkedHashMap<String, String> properties) {
        try {
            return dxApiProcessManager.propertiesPost(id, properties).toString();
        } catch (PropertiesException e) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(e.toString()).toString();
        }
    }

    @RequestMapping(value = "/sample/{id}/properties", method = RequestMethod.GET)
    public String samplePropertiesGet(
            @PathVariable Integer id) {
        try {
            return dxApiProcessManager.propertiesGet(id).toString();
        } catch (PropertiesException e) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(e.toString()).toString();
        }
    }

    @RequestMapping(value = "/sample/{id}/properties", method = RequestMethod.PUT)
    public String samplePropertiesPut(
            @PathVariable Integer id,
            @RequestBody LinkedHashMap<String, String> properties) {
        try {
            return dxApiProcessManager.propertiesPut(id, properties).toString();
        } catch (PropertiesException e) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(e.toString()).toString();
        }
    }

    @RequestMapping(value = "/sample/{id}/properties", method = RequestMethod.DELETE)
    public String samplePropertiesDelete(
            @PathVariable Integer id,
            @RequestBody LinkedHashMap<String, String> properties) {
        try {
            return dxApiProcessManager.propertiesDelete(id, properties).toString();
        } catch (PropertiesException e) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(e.toString()).toString();
        }
    }

}
