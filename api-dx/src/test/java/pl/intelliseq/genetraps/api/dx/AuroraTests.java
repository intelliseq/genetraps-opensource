package pl.intelliseq.genetraps.api.dx;

import lombok.extern.log4j.Log4j2;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.boot.test.context.SpringBootTest.WebEnvironment;
import org.springframework.test.context.junit4.SpringRunner;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

import java.sql.SQLException;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static pl.intelliseq.genetraps.api.dx.TestUser.ADMIN;
import static pl.intelliseq.genetraps.api.dx.TestUser.DEVIL;
import static pl.intelliseq.genetraps.api.dx.TestUser.PSYDUCK;


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
        Roles testRole = auroraDBManager.getUserPrivilegesToSample(1, 1);
        assertEquals(testRole, Roles.ADMIN);

        Map<String, Object> psyduck = auroraDBManager.getUserDetails(PSYDUCK.getId());
        assertEquals(psyduck.get("Username"), PSYDUCK.getUsername());
        assertEquals(psyduck.get("UserID"), PSYDUCK.getId());
        assertEquals(psyduck.get("Root"), 1);


        Map<String, Object> admin = auroraDBManager.getUserDetails(ADMIN.getUsername());
        assertEquals(admin.get("Username"), ADMIN.getUsername());
        assertEquals(admin.get("UserID"), ADMIN.getId());
        assertEquals(admin.get("Root"), 1);


        Map<String, Object> devil = auroraDBManager.getUserDetails(DEVIL.getUsername());
        assertEquals(devil.get("Username"), DEVIL.getUsername());
        assertEquals(devil.get("UserID"), DEVIL.getId());
        assertEquals(devil.get("Root"), 0);


        Map<Integer, Roles> adminRoles = auroraDBManager.getUserPrivileges(auroraDBManager.createNewSimpleUser("admin"));
        assertEquals(adminRoles.get(1), Roles.ADMIN);

    }

//    @Test
//    public void privileges(){
//        auroraDBManager.setUserPrivilegesToSample(PSYDUCK.toString().toLowerCase(), 5, Roles.ADMIN);
//    }

}
