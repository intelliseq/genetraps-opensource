package pl.intelliseq.genetraps.api.dx.helpers;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import lombok.extern.log4j.Log4j2;
import org.json.JSONObject;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.SqlParameter;
import org.springframework.jdbc.core.simple.SimpleJdbcCall;
import pl.intelliseq.genetraps.api.dx.Roles;
import pl.intelliseq.genetraps.api.dx.User;

import javax.annotation.PostConstruct;
import javax.sql.DataSource;
import java.sql.Types;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by intelliseq on 26/04/2017.
 */
@Log4j2
public class AuroraDBManager {

    @Autowired
    private Environment env;

    @Autowired
    private DataSource dataSource;

    @Autowired
    private FilesManager filesManager;

    private JdbcTemplate jdbcTemplate;

    @PostConstruct
    private void postConstruct() {
        jdbcTemplate = new JdbcTemplate(dataSource);
    }

    public List<String> getJobs(boolean withStatusNotYetSucceeded) {
        String sql;
        if(withStatusNotYetSucceeded)
            sql = "SELECT JobID FROM Jobs WHERE JobStatus = 0";
        else
            sql = "SELECT JobID FROM Jobs";
        // columns like: 1, 2, 3, ...
        // so not so much programmer..
        return jdbcTemplate.queryForList(sql, String.class);
    }

    public List<String> getJobsWithSampleID(String sampleId) {
        //TODO: add quotation marks around %s when sampleid changed to hash (string)
        String sql = String.format("SELECT JobID FROM Jobs WHERE SampleID = %s", sampleId);
        return jdbcTemplate.queryForList(sql, String.class);
    }

    public Integer checkWdlId(String wdlName) {
        String sql = String.format("SELECT WdlID FROM Wdls WHERE WdlName = \"%s\"", wdlName);
        return jdbcTemplate.query(sql, rs -> rs.next() ? rs.getInt("WdlID") : null);
    }

    public void putWdlToDB(String wdlName) {
        jdbcTemplate.update("INSERT INTO Wdls VALUES (?, ?)", 0, wdlName);
    }

    public void putJobToDB(String jobId, Integer userId, Integer wdlId, Integer jobStatus, Integer sampleId, JSONObject output) {
        jdbcTemplate.update("INSERT INTO Jobs VALUES (?, ?, ?, ?, ?, ?)", jobId, userId, wdlId, jobStatus, sampleId, output.toString());
    }

    /*
        User management
     */

    public User putUserToDB(User user, String password) {
        jdbcTemplate.update("INSERT INTO Users VALUES (?, ?, ?, ?)", 0, user.getLastName(), user.getFirstName(), user.getEmail());
        Integer userId = jdbcTemplate.queryForObject(String.format("SELECT UserID FROM Users WHERE Email = \"%s\"", user.getEmail()), Integer.class);
        user.setId(userId);
        jdbcTemplate.update("INSERT INTO Security VALUES (?, ?, ?, ?)", userId, user.getUserName(), password, user.getRoot() ? 1 : 0);
        return user;
    }

    public User getUserDetails(Integer userId) {
        SimpleJdbcCall jdbcCall = new SimpleJdbcCall(jdbcTemplate)
                .withProcedureName("GetUserDetailsByUserID")
                .declareParameters(new SqlParameter("UserID", Types.SMALLINT));

        jdbcCall.compile();

        Map<String, Object> inParamMap = new HashMap<>();
        inParamMap.put("UserID", userId);

        Map<String, Object> map = ((List<Map<String, Object>>) jdbcCall.execute(inParamMap).get("#result-set-1")).get(0);

        log.debug("User" + userId + " " + map.toString());
        log.debug(map.get("Root").getClass());

        return User.builder()
                .withId(userId)
                .withEmail((String) map.get("Email"))
                .withFirstName((String) map.get("FirstName"))
                .withLastName((String) map.get("LastName"))
                .withUserName((String) map.get("UserName"))
                .withRoot(map.get("Root").equals(1)).build();


    }

    public void getUsers() {
        jdbcTemplate.query("SELECT * FROM Users", (rs, rn) -> System.out.printf("%s:\t%s %s\n", rn, rs.getString("FirstName"), rs.getString("LastName")));
    }

    public Integer addUserToGroup(Integer userId, Integer groupId) {
        return jdbcTemplate.update("INSERT INTO UserGroups VALUES (?, ?)", userId, groupId);
    }

    /*
        Group managment
     */

    public JsonNode getUsersGroups(Integer userId) {
        return getUsersGroupsFromQuery(String.format(
                "SELECT G.* FROM Users AS U\n" +
                        "LEFT JOIN UserGroups AS UG\n" +
                        "ON U.UserID = UG.UserID\n" +
                        "LEFT JOIN Groups AS G\n" +
                        "ON UG.GroupID = G.GroupID\n" +
                        "WHERE U.UserID = %d;", userId
        ));
    }

    private JsonNode getUsersGroupsFromQuery(String query) {
        ArrayNode node = new ObjectMapper().createArrayNode();

        jdbcTemplate.query(query, (resultSet, rowNum) -> node.add(new ObjectMapper().createObjectNode()
                .put("GroupID", resultSet.getInt("GroupID"))
                .put("GroupName", resultSet.getString("GroupName"))
                .put("Root", resultSet.getInt("Root") == 1)));

        log.debug(node);

        return node;
    }

    public Integer createGroup(String groupName, Boolean root) {
        return jdbcTemplate.update("INSERT INTO Groups VALUES (?, ?, ?)", 0, groupName, root ? 1 : 0);
    }

    /*
        Priviledges
     */

    public Roles getUserPrivilegesToSample(Integer userID, Integer sampleID) {
        User user = getUserDetails(userID);

        if (user.getRoot()) {
            return Roles.ADMIN;
        }

        SimpleJdbcCall jdbcCall = new SimpleJdbcCall(jdbcTemplate)
                .withProcedureName("GetUserPrivilegesToSample")
                .declareParameters(new SqlParameter("UserID", Types.SMALLINT))
                .declareParameters(new SqlParameter("SampleID", Types.SMALLINT));

        jdbcCall.compile();

        Map<String, Object> inParamMap = new HashMap<>();
        inParamMap.put("UserID", user.getId());
        inParamMap.put("SampleID", sampleID);

        log.debug(jdbcCall.getInParameterNames());

        return Roles.valueOf(((List<Map<String, String>>) jdbcCall.execute(inParamMap).get("#result-set-1")).get(0).get("RoleName"));
    }

    private Map<Integer, Roles> getRootPrivileges() {
        return filesManager.getNumericDirectories().stream().collect(Collectors.toMap(e -> e, e -> Roles.ADMIN));
    }

    public Map<Integer, Roles> getUserPrivileges(Integer userID) {
        User user = getUserDetails(userID);

        if (user.getRoot()) {
            return getRootPrivileges();
        }
        SimpleJdbcCall jdbcCall = new SimpleJdbcCall(jdbcTemplate)
                .withProcedureName("GetUserPrivileges")
                .declareParameters(new SqlParameter("UserID", Types.SMALLINT));

        jdbcCall.compile();

        Map<String, Object> inParamMap = new HashMap<>();
        inParamMap.put("UserID", user.getId());

        log.debug(jdbcCall.getInParameterNames());


        //Changing something very bad to something less bad.
        Map<Integer, Roles> result = ((List<Map<String, Object>>) jdbcCall.execute(inParamMap).get("#result-set-1"))
                .stream()
                .collect(
                        Collectors.toMap(
                                e -> (Integer) e.get("SampleID"),
                                e -> Roles.valueOf(
                                        (String) e.get("RoleName")
                                )
                        )
                );

        log.info(result);

        return result;

    }


    public Integer setUserPrivilegesToSample(Integer userId, Integer sampleId, Roles role) {
        return jdbcTemplate.update("INSERT INTO UserPrivileges VALUES (?, ?, ?)", userId, sampleId, role.getId());
    }


}
