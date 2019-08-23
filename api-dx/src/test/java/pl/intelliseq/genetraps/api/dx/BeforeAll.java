package pl.intelliseq.genetraps.api.dx;

import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.testng.AbstractTestNGSpringContextTests;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.BeforeTest;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

import static pl.intelliseq.genetraps.api.dx.UserTest.PSYDUCK;

//@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = SpringBootTest.WebEnvironment.RANDOM_PORT)
@Log4j2
public class BeforeAll extends AbstractTestNGSpringContextTests {


    @Autowired
    private AuroraDBManager auroraDBManager;

    @BeforeTest
    public void beforeTest() throws Exception {
        super.springTestContextPrepareTestInstance();
        System.out.println("beforeTest");
        User psyduck = auroraDBManager.getUserDetails(PSYDUCK.getId());
        log.info(psyduck);

    }

    @BeforeMethod
    public void beforeMethod()
    {
        System.out.println("beforeMethod");
    }
}
