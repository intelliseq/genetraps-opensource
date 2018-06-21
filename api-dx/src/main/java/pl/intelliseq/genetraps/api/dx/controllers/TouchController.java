package pl.intelliseq.genetraps.api.dx.controllers;

import com.dnanexus.DXJob;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.ResponseBody;
import org.springframework.web.bind.annotation.RestController;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;

@RestController
public class TouchController {

    private Logger log = Logger.getLogger(TouchController.class);

    @Autowired
    private FilesManager filesManager;

    @Autowired
    private DxApiProcessManager processManager;

    @RequestMapping(value = "/touch", method = RequestMethod.GET)
    public String vep() {
        DXJob result = processManager.runTouch("-iname=test");
        log.info(result);
        return result.toString();
    }

    @RequestMapping(value = "/mkdir", method = RequestMethod.GET)
    @ResponseBody
    public String mkDir() {
        return String.format("{\"response\":%s}", filesManager.mkdir());
    }

}
