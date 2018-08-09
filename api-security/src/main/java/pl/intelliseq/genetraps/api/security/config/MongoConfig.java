package pl.intelliseq.genetraps.api.security.config;

import com.mongodb.MongoClient;
import com.mongodb.MongoClientURI;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.PropertySource;
import org.springframework.core.env.Environment;
import org.springframework.data.mongodb.config.AbstractMongoConfiguration;
import org.springframework.data.mongodb.repository.config.EnableMongoRepositories;


@Configuration
@EnableMongoRepositories()
@PropertySource(value={"classpath:application.properties"})
public class MongoConfig extends AbstractMongoConfiguration {
    @Value("${spring.data.mongodb.host}")
    private String ip;

    @Value("${spring.data.mongodb.database}")
    private String database;

    @Value("${spring.data.mongodb.authenticationDatabase}")
    private String authDb;

    @Autowired
    private Environment env;

    @Override
    protected String getDatabaseName() {
        return database;
    }

    @Override
    @Bean(name = "mongoClient")
    public MongoClient mongoClient(){
        String connection = String.format(
                "mongodb://%s:%s@%s/%s",
                env.getProperty("MONGO_ADMIN_NAME"),
                env.getProperty("MONGO_ADMIN_PASSWORD"),
                ip,
                authDb);
        return new MongoClient(new MongoClientURI(connection));
    }
}