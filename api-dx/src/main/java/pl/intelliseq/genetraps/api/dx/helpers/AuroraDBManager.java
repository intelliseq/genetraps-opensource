package pl.intelliseq.genetraps.api.dx.helpers;

import com.dnanexus.DXFile;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.module.SimpleModule;
import com.fasterxml.jackson.databind.node.ObjectNode;
import com.mysql.jdbc.DatabaseMetaData;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.SqlParameter;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
import org.springframework.jdbc.core.namedparam.SqlParameterSource;
import org.springframework.jdbc.core.simple.SimpleJdbcCall;

import javax.annotation.PostConstruct;
import javax.sql.DataSource;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by intelliseq on 26/04/2017.
 */
@Log4j2
public class AuroraDBManager {


    @Autowired
    private Environment env;

    @Autowired
    private DataSource dataSource;

    private JdbcTemplate jdbcTemplate;

    @PostConstruct
    private void postConstruct(){
        jdbcTemplate = new JdbcTemplate(dataSource);
    }

    public Map<String, Object> getUserPriviligeToSample(Integer userID, Integer sampleID) throws SQLException {
        SimpleJdbcCall jdbcCall = new SimpleJdbcCall(jdbcTemplate)
                .withProcedureName("GetUserPriviligeToSample")
                .declareParameters(new SqlParameter("UserID", Types.SMALLINT))
                .declareParameters(new SqlParameter("SampleID", Types.SMALLINT));

        jdbcCall.compile();

        Map<String, Object> inParamMap = new HashMap<>();
        inParamMap.put("UserID", userID);
        inParamMap.put("SampleID", sampleID);
        SqlParameterSource in = new MapSqlParameterSource(inParamMap);

        log.debug(jdbcCall.getInParameterNames());

        return(jdbcCall.execute(inParamMap));

    }

    public JsonNode getUserSimpleDetails(String username){

        String query = String.format("SELECT U.*, S.Username FROM " +
                "Security AS S " +
                "JOIN Users AS U ON S.UserID=U.USerID " +
                "WHERE S.Username = \"%s\";", username);

//        log.debug(query);
//
//        jdbcTemplate.query(query, (rs, rowNum) -> System.out.printf("Email: %s\n", rs.getString("Email")));
//
//        SimpleModule module = new SimpleModule();
//        module.addSerializer(new ResultSetSerializer());
//
//        ObjectMapper objectMapper = new ObjectMapper();
//        objectMapper.registerModule(module);

        ObjectNode node = jdbcTemplate.query(query, (resultSet, rowNum) -> new ObjectMapper().createObjectNode()
                .put("UserID", resultSet.getInt("UserID"))
                .put("LastName", resultSet.getString("LastName"))
                .put("FirstName", resultSet.getString("FirstName"))
                .put("Email", resultSet.getString("Email"))
                .put("Username", resultSet.getString("Username"))).get(0);

        log.debug(node);

        return node;

//        ResultSet resultSet = jdbcTemplate.query(query, (rs, rowNum) -> rs).get(0);
//
//        try {
//            log.debug(resultSet.getString("Email"));
//        } catch (SQLException e) {
//            throw new RuntimeException(e);
//        }
//
//        //TODO: Można napisać konwerter (https://stackoverflow.com/a/8120442)
//        JsonNode user = null;
//        try {
//            user = new ObjectMapper().createObjectNode()
//                    .put("UserID", resultSet.getInt("UserID"))
//                    .put("LastName", resultSet.getString("LastName"))
//                    .put("FirstName", resultSet.getString("FirstName"))
//                    .put("Email", resultSet.getString("Email"))
//                    .put("Username", username);
//        } catch (SQLException e) {
//            throw new RuntimeException(e);
//        }
//
//        return user;
    }

    public void getUsers(){
        jdbcTemplate.query("SELECT * FROM Users", (rs, rn) -> System.out.printf("%s:\t%s %s\n", rn, rs.getString("FirstName"), rs.getString("LastName")));
    }
}