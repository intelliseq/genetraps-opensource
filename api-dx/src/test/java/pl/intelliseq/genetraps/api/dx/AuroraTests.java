package pl.intelliseq.genetraps.api.dx;

import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.test.context.testng.AbstractTestNGSpringContextTests;
import org.testng.annotations.Test;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

import java.sql.SQLException;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static pl.intelliseq.genetraps.api.dx.TestUser.ADMIN;
import static pl.intelliseq.genetraps.api.dx.TestUser.DEVIL;
import static pl.intelliseq.genetraps.api.dx.TestUser.PSYDUCK;

//@RunWith(SpringRunner.class)
@SpringBootTest(webEnvironment = WebEnvironment.RANDOM_PORT)
@Log4j2
//@ActiveProfiles("test")
public class AuroraTests extends AbstractTestNGSpringContextTests {

    @Autowired
    private AuroraDBManager auroraDBManager;

    @Test
    public void auroraTest() throws SQLException {
        auroraDBManager.getUsers();
        Roles testRole = auroraDBManager.getUserPrivilegesToSample(1, 1);
        assertEquals(testRole, Roles.ADMIN);

        User psyduck = auroraDBManager.getUserDetails(PSYDUCK.getId());
        assertEquals(psyduck.getUserName(), PSYDUCK.getUsername());
        assertEquals(psyduck.getId(), PSYDUCK.getId());
        assertEquals(psyduck.getRoot(), true);


        User admin = auroraDBManager.getUserDetails(ADMIN.getId());
        assertEquals(admin.getUserName(), ADMIN.getUsername());
        assertEquals(admin.getId(), ADMIN.getId());
        assertEquals(admin.getRoot(), true);


        User devil = auroraDBManager.getUserDetails(DEVIL.getId());
        assertEquals(devil.getUserName(), DEVIL.getUsername());
        assertEquals(devil.getId(), DEVIL.getId());
        assertEquals(devil.getRoot(), false);


        Map<Integer, Roles> adminRoles = auroraDBManager.getUserPrivileges(ADMIN.getId());
        assertEquals(adminRoles.get(1), Roles.ADMIN);

    }

    @Test
    public void createUser(){
        User user = User.builder().withUserName("Test").withFirstName("Test").withLastName("Test").withEmail("test"+System.currentTimeMillis()+"@test.test").build();
        auroraDBManager.putUserToDB(user, "$2a$04$rIGVUoinjuyNn1TFh9FSA.YVLjfV7DfV52iabUUBR00QIIqQrddSe");
    }

//    @Test
//    public void privileges(){
//        auroraDBManager.setUserPrivilegesToSample(PSYDUCK.toString().toLowerCase(), 5, Roles.ADMIN);
//    }

}
