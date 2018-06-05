package pl.intelliseq.genetraps.api.dx;

import com.dnanexus.DXContainer;
import org.apache.log4j.Logger;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.core.env.Environment;
import org.springframework.test.context.junit4.SpringRunner;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;

import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.assertj.core.api.Assertions.assertThat;
import static org.junit.Assert.assertEquals;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
public class ApiDxApplicationTests {

	Logger log = Logger.getLogger(ApiDxApplicationTests.class);

	@Autowired
	private TestRestTemplate restTemplate;

	@Autowired
	private FilesManager filesManager;

	@Autowired
	private Environment env;

	@Test
	public void contextLoads() {
	}

	@Test
	public void mkdirTest() throws InterruptedException {
		class MkDirRunnable implements Runnable{
			private int id;

			public MkDirRunnable(int id){
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

		IntStream.range(0,10).forEach(i -> threadPoolExecutor.execute(new MkDirRunnable(i)));

		while (threadPoolExecutor.getTaskCount()!=threadPoolExecutor.getCompletedTaskCount()){
			log.debug("count="+threadPoolExecutor.getTaskCount()+","+threadPoolExecutor.getCompletedTaskCount());
			Thread.sleep(5000);
		}
		threadPoolExecutor.shutdown();
		threadPoolExecutor.awaitTermination(60, TimeUnit.SECONDS);

        assertEquals(filesManager.getNumericDirectories().size() - size, 10);

	}
	
}
