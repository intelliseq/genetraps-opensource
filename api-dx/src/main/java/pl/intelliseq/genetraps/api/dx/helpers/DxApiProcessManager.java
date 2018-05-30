package pl.intelliseq.genetraps.api.dx.helpers;

import com.dnanexus.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import pl.intelliseq.genetraps.api.dx.IProcessManager;

import java.util.HashMap;

public class DxApiProcessManager {

    @Autowired
    Environment env;

    public DXJob runUrlFetch(String inputUrl, String sampleNumber, String... tags) {
        var input = new HashMap<>();
        input.put("url", inputUrl);
        input.put("tags", tags);

        return DXApp.getInstance("app-FF6b2Bj9pQGXV19347j6fJPX")
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setFolder(String.format("/samples/%s/rawdata", sampleNumber))
                .setInput(input)
                .run();
    }

    public void runMkDir(Integer sampleNumber){
        runMkDir(sampleNumber.toString());
    }

    public void runMkDir(String sampleNumber){
        DXContainer.getInstance(env.getProperty("dx-project")).newFolder("/samples/"+sampleNumber);
    }

//    public DXJob runTouch(String fileId){
//
//        return DXApplet.getInstance("applet-FG6V7Y0045kFGFvy5xf68X42")
//    }

    public DXJob runFastqc(String fileId){
        var input = new HashMap<>();
        input.put("fastq_file", DXFile.getInstance(fileId));

        return DXApplet.getInstance("applet-FG6V7qj045k7XVk43BZ34Fy8")
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setInput(input)
                .run();
    }
}
