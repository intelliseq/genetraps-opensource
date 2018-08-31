package pl.intelliseq.genetraps.api.dx.helpers;

import com.mysql.jdbc.DatabaseMetaData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
import org.springframework.jdbc.core.simple.SimpleJdbcCall;

import javax.annotation.PostConstruct;
import javax.sql.DataSource;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Map;

/**
 * Created by intelliseq on 26/04/2017.
 */
public class AuroraDBManager {


    @Autowired
    private Environment env;

    @Autowired
    private DataSource dataSource;

    private JdbcTemplate jdbcTemplate;

    private Logger logger = LoggerFactory.getLogger(AuroraDBManager.class);

    @PostConstruct
    private void postConstruct(){
        jdbcTemplate = new JdbcTemplate(dataSource);
    }

    public Map<String, Object> getUserPriviligeToSample(Integer userID, Integer sampleID) throws SQLException {



        ResultSet rs = dataSource.getConnection().getMetaData().getProcedures("genetraps_security", "genetraps_security", "%");

        while (rs.next()) {
            String spName = rs.getString("PROCEDURE_NAME");
            int spType = rs.getInt("PROCEDURE_TYPE");
            System.out.println("Stored Procedure Name: " + spName);
            if (spType == DatabaseMetaData.procedureReturnsResult) {
                System.out.println("procedure Returns Result");
            } else if (spType == DatabaseMetaData.procedureNoResult) {
                System.out.println("procedure No Result");
            } else {
                System.out.println("procedure Result unknown");
            }

        }

        return null;
//
//        SimpleJdbcCall jdbcCall = new SimpleJdbcCall(dataSource)
//                .withCatalogName("genetraps_security")
//                .withSchemaName("genetraps_security")
//                .withProcedureName("GetUserPriviligeToSample");
//
//        return(jdbcCall.execute(new MapSqlParameterSource()
//                .addValue("UserID", userID)
//                .addValue("SampleID", sampleID)));

    }

    public void getUsers(){
        jdbcTemplate.query("SELECT * FROM Users", (rs, rn) -> System.out.printf("%s:\t%s %s\n", rn, rs.getString("FirstName"), rs.getString("LastName")));
    }
}
