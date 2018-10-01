package pl.intelliseq.genetraps.api.dx;

import com.dnanexus.DXContainer;
import com.dnanexus.DXEnvironment;
import com.dnanexus.DXHTTPRequest;
import com.dnanexus.DXJSON;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.log4j.Log4j2;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.core.env.Environment;
import org.springframework.test.context.ActiveProfiles;
import org.springframework.test.context.junit4.SpringRunner;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;

import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import static org.junit.Assert.assertEquals;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
@ActiveProfiles("test")
public class ApiDxTest {

    @Autowired
    FilesManager filesManager;


    @Autowired
    Environment env;

    @Autowired
    DxApiProcessManager dxApiProcessManager;

    @Test
    public void singleMkdirTest() {
        Integer folder = filesManager.getLowestFreeIndex();
        dxApiProcessManager.runMkDir(folder);
    }

    //    @Test
    public void descriptiontesting() {
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

    @Test
    public void manyMkdirTest() throws InterruptedException {
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

        DXContainer.getInstance(env.getProperty("dx-project")).removeFolder("/samples/" + folderToDelete);

        assertEquals(filesManager.mkdir(), folderToDelete);

    }

}
