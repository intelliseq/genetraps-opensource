package pl.intelliseq.genetraps.api.dx;

import org.apache.log4j.Logger;
import org.hamcrest.Matchers;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.MediaType;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.util.LinkedMultiValueMap;
import org.springframework.util.MultiValueMap;
import pl.intelliseq.genetraps.api.dx.enums.JobStates;
import pl.intelliseq.genetraps.api.dx.exceptions.IseqParseException;
import pl.intelliseq.genetraps.api.dx.models.IseqJSON;

import static org.junit.Assert.assertThat;
import static org.junit.Assert.fail;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
public class UploadTests {

	Logger log = Logger.getLogger(UploadTests.class);

	@Autowired
	private TestRestTemplate restTemplate;

    public String sampleLeft = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.1.fq.gz";
    public String sampleRight = "http://resources.intelliseq.pl/kamilant/test-data/fastq/capn3.2.fq.gz";

    public IseqJSON waitUntilJobIsDone(String jobId){
        String response;
        try{
            JobStates jobState;
            IseqJSON iseqJSON;
            do {
                response = this.restTemplate.getForObject("/describe/{id}", String.class, jobId);
                iseqJSON = new IseqJSON(response);
                String state = iseqJSON.getString("state");
                jobState = JobStates.getEnum(state);
                Thread.sleep(5000);
            }while (jobState != JobStates.DONE);
            return iseqJSON;
        }catch (InterruptedException e){
            throw new RuntimeException(e);
        }
    }

    public IseqJSON upload(String sampleUrl, Integer sampleNumber, String... tags){
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
        uploadValueMap.add("url", sampleUrl);
        uploadValueMap.add("sampleNumber", sampleNumber.toString());
        for(String tag:tags){
            uploadValueMap.add("tags", tag);
        }

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);

        String response = this.restTemplate.postForObject("/upload", uploadEntity, String.class);
        String jobId = new IseqJSON(response).getString("id");

        assertThat(jobId, Matchers.matchesPattern("job-\\w*"));

        IseqJSON jobResponse = waitUntilJobIsDone(jobId);

        String fileId = jobResponse.getIseqJSON("output").getIseqJSON("file").getString("$dnanexus_link");
        return new IseqJSON(this.restTemplate.getForObject("/describe/{id}", String.class, fileId));
    }

    public IseqJSON fastqc(String fileId){
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
        uploadValueMap.add("fileId", fileId);

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);

        String response = this.restTemplate.postForObject("/fastqc", uploadEntity, String.class);
        String jobId = new IseqJSON(response).getString("id");

        assertThat(jobId, Matchers.matchesPattern("job-\\w*"));

        return waitUntilJobIsDone(jobId);
    }

    public IseqJSON bwa(String left, String right, String outputFolder){
        HttpHeaders uploadHeaders = new HttpHeaders();
        uploadHeaders.setContentType(MediaType.APPLICATION_FORM_URLENCODED);

        MultiValueMap<String, String> uploadValueMap = new LinkedMultiValueMap<String, String>();
        uploadValueMap.add("left", left);
        uploadValueMap.add("right", right);
        uploadValueMap.add("outputFolder", outputFolder);

        HttpEntity<MultiValueMap<String, String>> uploadEntity = new HttpEntity<MultiValueMap<String, String>>(uploadValueMap, uploadHeaders);

        String response = this.restTemplate.postForObject("/bwa", uploadEntity, String.class);
        String jobId = new IseqJSON(response).getString("id");

        assertThat(jobId, Matchers.matchesPattern("job-\\w*"));

        return waitUntilJobIsDone(jobId);
    }

    @Test
    public void upload(){
        IseqJSON file1 = upload(sampleLeft, 1, "left");
        String file1Id = file1.getString("id");
        String file1Folder = file1.getString("folder").replace("/rawdata", "");
        String file2Id = upload(sampleRight, 1, "right").getString("id");

        log.info(bwa(file1Id, file2Id, file1Folder));

        log.info(fastqc(file1Id));
        log.info(fastqc(file2Id));
    }


	
}
