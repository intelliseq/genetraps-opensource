package pl.intelliseq.genetraps.api.dx;

import com.dnanexus.DXContainer;
import com.dnanexus.DXFile;
import com.fasterxml.jackson.databind.JsonNode;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.core.env.Environment;
import org.springframework.core.io.ClassPathResource;
import org.springframework.http.HttpEntity;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpMethod;
import org.springframework.http.ResponseEntity;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.test.context.testng.AbstractTestNGSpringContextTests;
import org.testng.annotations.Test;

import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.text.IsEqualIgnoringCase.equalToIgnoringCase;
import static org.junit.Assert.*;
import static pl.intelliseq.genetraps.api.dx.TestUser.PSYDUCK;


//@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
//@ActiveProfiles("test")
public class FileTests extends AbstractTestNGSpringContextTests {

    @Autowired
    Environment env;

    @Test
    public void main() throws IOException {

        DXFile.newFile()
                .setProject(DXContainer.getInstance(env.getProperty("dx-project")))
                .upload(
                new FileInputStream(new ClassPathResource("uploadTestFile.txt").getFile())
        ).build();
    }
}