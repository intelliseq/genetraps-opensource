package pl.intelliseq.genetraps.api.security;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.core.env.Environment;
import org.springframework.jdbc.datasource.DriverManagerDataSource;

import javax.sql.DataSource;

@Configuration
public class BeansConfiguration {

    @Autowired
    private Environment env;

    @Bean
    public DataSource dataSource() {
        DriverManagerDataSource ds = new DriverManagerDataSource();
        ds.setDriverClassName(com.mysql.jdbc.Driver.class.getName());
        ds.setUrl("jdbc:mysql://genetraps-provisioned.cluster-cvdvrgz3bmic.us-east-1.rds.amazonaws.com:3306/genetraps_security");
        ds.setUsername("genetraps-client");
        ds.setPassword(env.getProperty("AURORA_GENETRAPS_CLIENT_PASSWD"));
        return ds;
    }

    @Bean
    AuroraDBManager auroraDBManager(){
        return new AuroraDBManager();
    }
}
