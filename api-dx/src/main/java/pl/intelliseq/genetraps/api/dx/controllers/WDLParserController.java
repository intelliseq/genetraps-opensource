package pl.intelliseq.genetraps.api.dx.controllers;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;
import pl.intelliseq.genetraps.api.dx.helpers.WDLParserManager;

@RestController
public class WDLParserController {

    @Autowired
    WDLParserManager WDLParserManager;

    @RequestMapping(value = "/wdl/{name}/info", method = RequestMethod.GET)
    public String GetSpecificWDL(
            @PathVariable String name) {
        return WDLParserManager.GetData(name).toString();
    }

    @RequestMapping(value = "/wdl/info", method = RequestMethod.GET)
    public String GetAllWDLs() {
        return WDLParserManager.GetData().toString();
    }

    @RequestMapping(value = "/wdl/update", method = RequestMethod.GET)
    public void UpdateWDLData() {
        WDLParserManager.CollectData();
    }
}
