package pl.intelliseq.genetraps.api.dx;

import org.apache.log4j.Logger;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.test.context.junit4.SpringRunner;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;
import pl.intelliseq.genetraps.api.dx.helpers.ProcessManager;

import java.io.IOException;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
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
	private ProcessManager processManager;

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

	@Test
	public void retrieveLastDirIndex(){
        Integer lowest = filesManager.getLowestFreeIndex();
        if(lowest > 1){
            Integer selected = lowest/2;
            filesManager.resetCounter();
            log.info(String.format("Lowest: %d Selected: %d", lowest, selected));
            processManager.runCommand("dx rmdir "+selected);
            log.info("Check if lowest index is "+selected);
            assertEquals(filesManager.getLowestFreeIndex(), selected);
            assertEquals(filesManager.mkdir(), selected);
            log.info("Check if lowest index is back "+lowest);
            assertEquals(filesManager.getLowestFreeIndex(), lowest);
        }

	}

	@Test
	public void testVep() {
		
		//MultiValueMap<String, Object> map = new LinkedMultiValueMap<String, Object>();
		//map.add("variant", "rs669");
		
		String body = this.restTemplate.getForObject("/touch", String.class);
		//String body = this.restTemplate.postForObject("/touch", map, String.class);
		System.out.println(body);
		assertThat(body.contains("Job ID"));
	}
	
//	//@Test
//	public void dxTest() throws IOException, InterruptedException {
//    	String result = processManager.runCommand("printf \"Y\\n\" | dx run touch -iname=test");
//
//    	//log.info("Result: " + result);
//    	//log.info("Result: " + this.getJobId(result + "\ntadam"));
//
//    	String jobId = this.getJobId(result);
//
//    	//log.info(DxJob.getDxJobById(jobId));
//
//    	for (int i = 1; i <= 30; i++) {
//    		log.info("i: " + i);
//    		Thread.sleep(500);
//    		log.info(DxJob.getDxJobById(jobId));
//    	}
//	}
//
//	private String getJobId(String result) {
//		try {
//			int indexOfJobId = result.indexOf("Job ID") + 8;
//			int endIndexOfJobId = result.substring(indexOfJobId).indexOf("\n");
//			if (endIndexOfJobId == -1) {
//				return result.substring(indexOfJobId);
//			}
//			return result.substring(indexOfJobId, indexOfJobId + endIndexOfJobId);
//		} catch (Exception e) {
//			return null;
//		}
//	}
	
}
