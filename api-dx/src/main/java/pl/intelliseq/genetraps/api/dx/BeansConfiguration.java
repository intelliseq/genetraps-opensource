package pl.intelliseq.genetraps.api.dx;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;

/**
 * Created by intelliseq on 07/12/2017.
 */
@Configuration
public class BeansConfiguration {
    @Bean
    FilesManager filesManager(){
        return new FilesManager();
    }
}
