package pl.intelliseq.genetraps.api.security;

import com.mysql.jdbc.DatabaseMetaData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.RowMapper;
import org.springframework.security.core.authority.SimpleGrantedAuthority;
import org.springframework.security.core.userdetails.User;

import javax.annotation.PostConstruct;
import javax.sql.DataSource;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.Collections;
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

    public User getUser(String username){
        return jdbcTemplate.query(
                String.format("SELECT S.* FROM " +
                        "Security AS S " +
                        "JOIN Users AS U ON S.UserID=U.UserID " +
                        "WHERE S.Username = \"%s\" OR U.Email = \"%s\";", username, username),
                (rs, rowNum) -> new User(rs.getString("Username"), rs.getString("Password"), Collections.singletonList(new SimpleGrantedAuthority("user")))).get(0);
    }
}
