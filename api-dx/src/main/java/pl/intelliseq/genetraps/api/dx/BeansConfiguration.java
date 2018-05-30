package pl.intelliseq.genetraps.api.dx;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;
import pl.intelliseq.genetraps.api.dx.helpers.ProcessManager;

/**
 * Created by intelliseq on 07/12/2017.
 */
@Configuration
public class BeansConfiguration {
    @Bean
    FilesManager filesManager(){
        return new FilesManager();
    }
    @Bean
    ProcessManager processManager(){ return new ProcessManager(); }
    @Bean
    DxApiProcessManager dxApiProcessManager(){
        return new DxApiProcessManager();
    }
}
