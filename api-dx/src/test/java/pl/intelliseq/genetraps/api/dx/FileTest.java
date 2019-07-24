package pl.intelliseq.genetraps.api.dx;

import com.dnanexus.DXContainer;
import com.dnanexus.DXFile;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.core.env.Environment;
import org.springframework.core.io.ClassPathResource;
import org.springframework.test.context.testng.AbstractTestNGSpringContextTests;
import org.testng.annotations.Test;

import java.io.FileInputStream;
import java.io.IOException;

import static org.hamcrest.CoreMatchers.is;


//@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
//@ActiveProfiles("test")
public class FileTest extends AbstractTestNGSpringContextTests {

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