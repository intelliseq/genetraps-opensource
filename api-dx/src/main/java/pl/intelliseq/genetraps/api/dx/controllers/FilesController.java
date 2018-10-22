

import pl.intelliseq.genetraps.api.dx.exceptions.PropertiesException;
import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.security.oauth2.provider.OAuth2Authentication;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import pl.intelliseq.genetraps.api.dx.Roles;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;

import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
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
            @RequestParam String sampleId,
            @RequestParam String... tag) {
        log.info("upload");
        log.debug(Arrays.toString(tag));
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());
        log.debug(userId);

        return new ObjectMapper().createObjectNode().put("id", processManager.runUrlFetch(url, sampleId, tag).getId()).toString();
    }

    // document me change in master:readme
    @RequestMapping(value = "/uploadfile", method = RequestMethod.POST)
    public String uploadfile(
            @RequestParam(value = "file") MultipartFile file,
            @RequestParam(value = "sampleId") int sampleId,
            @RequestParam(value = "newfilename", required = false) String newfilename,
            @RequestParam(value = "tag", required = false) List<String> tags) {

        try {
            return new ObjectMapper().createObjectNode().put("id", processManager.runUploadFile(file, sampleId, newfilename, tags)).toString();
        } catch (IOException e) {
            return new ObjectMapper().createObjectNode().put("id", "").toString();
        }
    }

    @RequestMapping(value = "/mkdir", method = RequestMethod.GET)
    @ResponseBody
    public String mkDir() {
        return String.format("{\"response\":%s}", filesManager.mkdir());
    }

    // document me change in master:endpoints,readme
    @RequestMapping(value = "/sample/{id}/describe", method = RequestMethod.GET)
    public String describe(
            @PathVariable String id) {
        return processManager.JSONDescribe(id).toString();
    }

    // document me change in master:endpoints,readme
    @RequestMapping(value = "/sample/{id}/ls", method = RequestMethod.GET)
    public String sampleLs(@PathVariable("id") int sampleId,
                            @RequestParam(required = false, defaultValue = "false") boolean byNames) {
        return processManager.sampleLs(sampleId, byNames).toString();
    }

    @RequestMapping(value = "/sample/{id}/properties", method = RequestMethod.POST)
    public String samplePropertiesPost(@PathVariable("id") int sampleId,
                                @RequestBody LinkedHashMap<String, String> properties) {
        try {
            return processManager.propertiesPost(sampleId, properties).toString();
        } catch(PropertiesException e) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(e.toString()).toString();
        }
    }

    @RequestMapping(value = "/sample/{id}/properties", method = RequestMethod.GET)
    public String samplePropertiesGet(@PathVariable("id") int sampleId) {
        try {
            return processManager.propertiesGet(sampleId).toString();
        } catch (PropertiesException e) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(e.toString()).toString();
        }
    }

    @RequestMapping(value = "/sample/{id}/properties", method = RequestMethod.PUT)
    public String samplePropertiesPut(@PathVariable("id") int sampleId,
                                @RequestBody LinkedHashMap<String, String> properties) {
        try {
            return processManager.propertiesPut(sampleId, properties).toString();
        } catch (PropertiesException e) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(e.toString()).toString();
        }
    }

    @RequestMapping(value = "/sample/{id}/properties", method = RequestMethod.DELETE)
    public String samplePropertiesDelete(@PathVariable("id") int sampleId,
                                @RequestBody LinkedHashMap<String, String> properties) {
        try {
            return processManager.propertiesDelete(sampleId, properties).toString();
        } catch (PropertiesException e) {
            return ResponseEntity.status(HttpStatus.BAD_REQUEST).body(e.toString()).toString();
        }
    }

}
