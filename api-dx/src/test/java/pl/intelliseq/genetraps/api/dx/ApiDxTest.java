package pl.intelliseq.genetraps.api.dx;

import com.dnanexus.*;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.MappingJsonFactory;
import org.apache.log4j.Logger;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.core.env.Environment;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.MediaType;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.util.LinkedMultiValueMap;
import org.springframework.util.MultiValueMap;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;
import pl.intelliseq.genetraps.api.dx.helpers.ProcessManager;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
public class ApiDxTest {

    @Autowired
    FilesManager filesManager;

    @Autowired
    ProcessManager processManager;

    @Autowired
    Environment env;

    @Autowired
    DxApiProcessManager dxApiProcessManager;

	Logger log = Logger.getLogger(ApiDxTest.class);

//    @Test
    public void testDXAPICustomEnvironment() throws IOException {
        DXEnvironment env = DXEnvironment.Builder.fromDefaults().build();
        JsonNode input =
                (JsonNode) (new MappingJsonFactory().createParser("{}").readValueAsTree());
        JsonNode responseJson = DXAPI.systemFindDataObjects(input, JsonNode.class, env);
        Assert.assertEquals(responseJson.isObject(), true);
        DXEnvironment bogusContextEnv =
                DXEnvironment.Builder.fromDefaults().build();
        JsonNode output = DXAPI.systemFindDataObjects(input, JsonNode.class, bogusContextEnv);
        JsonNode output2 = DXAPI.projectDescribe("project-F5qXzZ8045k4x7V28ykjV2Gy", bogusContextEnv);
        System.out.println(output2);
        /*try {
            JsonNode output = DXAPI.systemFindDataObjects(input, JsonNode.class, bogusContextEnv);
            Assert.fail("Expected request with bogus token to throw InvalidAuthentication");
        } catch (InvalidAuthenticationException e) {
            // Expected
        }*/
    }

    @Test
    public void mkdirTest(){
        Integer folder = filesManager.getLowestFreeIndex();
        dxApiProcessManager.runMkDir(folder);
    }

//    @Test
    public void touch(){
        String sampleLeft = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.1.fq.gz";
        String[] tags = new String[]{"left", "one", "two"};

        var input = new HashMap<>();
        input.put("url", sampleLeft);
        input.put("tags", tags);

        DXApp.getInstance("app-FF6b2Bj9pQGXV19347j6fJPX")
                .newRun()
                .setProject(DXProject.getInstance(env.getProperty("dx-project")))
                .setFolder("/samples/7/rawdata")
                .setInput(input)
                .run();
    }

    @Test
    public void uploadTest(){
        String sampleLeft = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.1.fq.gz";
        String[] tags = new String[]{"left", "one", "two"};

        System.out.println(dxApiProcessManager.runUrlFetch(sampleLeft, "7", tags));
    }

    @Test
    public void fastqcTest(){
        dxApiProcessManager.runFastqc("file-FG2YvP00FjKBvJGx0x7BFzF3");
    }
	
}
