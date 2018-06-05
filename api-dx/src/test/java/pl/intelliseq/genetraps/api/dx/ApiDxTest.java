package pl.intelliseq.genetraps.api.dx;

import org.apache.log4j.Logger;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.core.env.Environment;
import org.springframework.test.context.junit4.SpringRunner;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;



@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
public class ApiDxTest {

    @Autowired
    FilesManager filesManager;


    @Autowired
    Environment env;

    @Autowired
    DxApiProcessManager dxApiProcessManager;

	Logger log = Logger.getLogger(ApiDxTest.class);

    @Test
    public void mkdirTest(){
        Integer folder = filesManager.getLowestFreeIndex();
        dxApiProcessManager.runMkDir(folder);
    }

//    @Test
//    public void touch(){
//        String sampleLeft = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.1.fq.gz";
//        String[] tags = new String[]{"left", "one", "two"};
//
//        var input = new HashMap<>();
//        input.put("url", sampleLeft);
//        input.put("tags", tags);
//
//        DXApp.getInstance("app-FF6b2Bj9pQGXV19347j6fJPX")
//                .newRun()
//                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
//                .setFolder("/samples/7/rawdata")
//                .setInput(input)
//                .run();
//    }

}
