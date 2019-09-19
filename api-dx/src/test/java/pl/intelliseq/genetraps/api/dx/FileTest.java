package pl.intelliseq.genetraps.api.dx;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.test.context.TestPropertySource;
import org.springframework.test.context.testng.AbstractTestNGSpringContextTests;
import org.testng.annotations.Test;
import pl.intelliseq.genetraps.api.dx.helpers.aws_manager.AWSApiProcessManagerImpl;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;

import static org.testng.Assert.assertEquals;


//@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@TestPropertySource("/application.properties")
@Log4j2
//@ActiveProfiles("test")
public class FileTest extends AbstractTestNGSpringContextTests {

    @Autowired
    private AWSApiProcessManagerImpl awsApiProcessManagerImpl;

    @Autowired
    private FilesManager filesManager;

    final private AmazonS3 s3Client = AmazonS3ClientBuilder.defaultClient();

    @Value("${bucket.default}")
    private String bucket;

    @Value("${samples.folder}")
    private String samplesFolder;

    @Value("${test.file.name}")
    private String testFileName;

//    @Test
//    public void main() throws IOException {
//
//        DXFile.newFile()
//                .setProject(DXContainer.getInstance(env.getProperty("dx-project")))
//                .upload(
//                new FileInputStream(new ClassPathResource("uploadTestFile.txt").getFile())
//        ).build();
//    }

    @Test
    public void testDeleteFile() throws InterruptedException {

        Integer id = filesManager.getLowestFreeIndex(filesManager.getNumericDirectories());
        Integer sampleId = awsApiProcessManagerImpl.runCreateSample(id);
        s3Client.putObject(bucket, String.format("%s/%s/test", samplesFolder, sampleId), "testDeleteFile");
        assertEquals(awsApiProcessManagerImpl.runSampleLs(sampleId, "").length(), 1);
        log.info(awsApiProcessManagerImpl.runDeleteFile(sampleId, String.format("%s/%s/%s", samplesFolder, sampleId, testFileName)));
        assertEquals(awsApiProcessManagerImpl.runSampleLs(sampleId, "").length(), 0);
        log.info(awsApiProcessManagerImpl.runDeleteFile(sampleId, String.format("%s/%s/", samplesFolder, sampleId)));
    }
}