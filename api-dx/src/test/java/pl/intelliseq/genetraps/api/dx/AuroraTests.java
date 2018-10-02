package pl.intelliseq.genetraps.api.dx;

import lombok.extern.log4j.Log4j2;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.boot.test.web.client.TestRestTemplate;
import org.springframework.test.context.ActiveProfiles;
import org.springframework.test.context.junit4.SpringRunner;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

import java.sql.SQLException;

import static org.junit.Assert.assertEquals;


@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
//@ActiveProfiles("test")
public class AuroraTests {

    @Autowired
    private AuroraDBManager auroraDBManager;

    @Test
    public void auroraTest() throws SQLException {
        auroraDBManager.getUsers();
        log.info(auroraDBManager.getUserPriviligeToSample(1, 1));
    }

}
