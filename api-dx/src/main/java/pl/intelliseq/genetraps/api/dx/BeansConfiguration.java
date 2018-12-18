package pl.intelliseq.genetraps.api.dx;

import com.google.common.base.Predicate;
import com.google.common.base.Predicates;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.core.env.Environment;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.springframework.jdbc.datasource.DriverManagerDataSource;
import org.springframework.security.oauth2.provider.token.TokenStore;
import org.springframework.security.oauth2.provider.token.store.JwtAccessTokenConverter;
import org.springframework.security.oauth2.provider.token.store.JwtTokenStore;
import org.springframework.util.FileCopyUtils;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;
import springfox.documentation.builders.PathSelectors;
import springfox.documentation.builders.RequestHandlerSelectors;
import springfox.documentation.spi.DocumentationType;
import springfox.documentation.spring.web.plugins.Docket;
import springfox.documentation.swagger2.annotations.EnableSwagger2;

import javax.sql.DataSource;
import java.io.IOException;

/**
 * Created by intelliseq on 07/12/2017.
 */
@Configuration
@EnableSwagger2
public class BeansConfiguration {
    @Bean
    FilesManager filesManager() {
        return new FilesManager();
    }

    @Bean
    DxApiProcessManager dxApiProcessManager() {
        return new DxApiProcessManager();
    }

    @Bean
    AuroraDBManager auroraDBManager() {
        return new AuroraDBManager();
    }

    // JWT Configuration
    @Autowired
    private JwtAccessTokenConverter jwtAccessTokenConverter;

    @Autowired
    private Environment env;

    @Bean
    @Qualifier("tokenStore")
    public TokenStore tokenStore() {
        return new JwtTokenStore(jwtAccessTokenConverter);
    }

    @Bean
    public DataSource dataSource() {
        DriverManagerDataSource ds = new DriverManagerDataSource();
        ds.setDriverClassName(com.mysql.jdbc.Driver.class.getName());
        ds.setUrl("jdbc:mysql://genetraps-provisioned.cluster-cvdvrgz3bmic.us-east-1.rds.amazonaws.com:3306/genetraps_security?noAccessToProcedureBodies=true");
        ds.setUsername("genetraps-client");
        ds.setPassword(env.getProperty("AURORA_GENETRAPS_CLIENT_PASSWD"));
        return ds;
    }

    @Bean
    protected JwtAccessTokenConverter jwtTokenEnhancer() {
        JwtAccessTokenConverter converter = new JwtAccessTokenConverter();
        Resource resource = new ClassPathResource("public.cert");
        String publicKey = null;
        try {
            publicKey = new String(FileCopyUtils.copyToByteArray(resource.getInputStream()));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        converter.setVerifierKey(publicKey);
        return converter;
    }

    @Bean
    public Docket api() {
        return new Docket(DocumentationType.SWAGGER_2)
                .select()
                .apis(RequestHandlerSelectors.any())
                .paths(paths())
                .build();
    }

    private Predicate<String> paths() {
        return Predicates.not(PathSelectors.regex("/(basic-)?error(-controller)?.*"));
    }
}
