package pl.intelliseq.genetraps.api.dx;

import com.dnanexus.DXContainer;
import com.dnanexus.DXEnvironment;
import com.dnanexus.DXHTTPRequest;
import com.dnanexus.DXJSON;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.core.env.Environment;
import org.springframework.test.context.ActiveProfiles;
import org.springframework.test.context.testng.AbstractTestNGSpringContextTests;
import org.testng.annotations.Test;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;

import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;


//@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
@ActiveProfiles("test")
public class ApiDxTest extends AbstractTestNGSpringContextTests {

    @Autowired
    FilesManager filesManager;


    @Autowired
    Environment env;

    @Autowired
    DxApiProcessManager dxApiProcessManager;

    @Test
    public void testSingleMkdir() {

        Integer folder = filesManager.getLowestFreeIndex(filesManager.getNumericDirectories());
//        dxApiProcessManager.runMkDir(folder);
        assertTrue(folder != null);
    }

    //    @Test
    public void testDescription() {
        var objectId = "job-FGP2j6Q0pqjbJZ7kKjJVJz8k";
        var out = DXJSON.safeTreeToValue(
                new DXHTTPRequest(DXEnvironment.create()).request("/" + objectId + "/" + "describe",
                        new ObjectMapper().createObjectNode(), DXHTTPRequest.RetryStrategy.SAFE_TO_RETRY), JsonNode.class);

        log.info(out);

        objectId = "file-FGbb8Q0021X1X4V2G1v3vv8G";
        out = DXJSON.safeTreeToValue(
                new DXHTTPRequest(DXEnvironment.create()).request("/" + objectId + "/" + "describe",
                        new ObjectMapper().createObjectNode(), DXHTTPRequest.RetryStrategy.SAFE_TO_RETRY), JsonNode.class);

        log.info(out);
    }

//    @Test
    public void testManyMkdir() throws InterruptedException {
        class MkDirRunnable implements Runnable {
            private int id;

            public MkDirRunnable(int id) {
                this.id = id;
            }

            @Override
            public void run() {
                log.info(id);
                filesManager.mkdir();
            }
        }

        Integer size = filesManager.getNumericDirectories().size();

        ThreadPoolExecutor threadPoolExecutor = (ThreadPoolExecutor) Executors.newFixedThreadPool(20);

        IntStream.range(0, 10).forEach(i -> threadPoolExecutor.execute(new MkDirRunnable(i)));

        while (threadPoolExecutor.getTaskCount() != threadPoolExecutor.getCompletedTaskCount()) {
            log.debug("count=" + threadPoolExecutor.getTaskCount() + "," + threadPoolExecutor.getCompletedTaskCount());
            Thread.sleep(5000);
        }
        threadPoolExecutor.shutdown();
        threadPoolExecutor.awaitTermination(60, TimeUnit.SECONDS);

        assertEquals(filesManager.getNumericDirectories().size() - size, 10);

        Integer folderToDelete = filesManager.getNumericDirectories().size() / 2;

        log.debug("Folder to delete: "+folderToDelete);

        DXContainer.getInstance(env.getProperty("dx-project")).removeFolder("/samples/" + folderToDelete, true);

        assertEquals(filesManager.mkdir(), folderToDelete);

    }

}
