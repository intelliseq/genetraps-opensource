package pl.intelliseq.genetraps.api.dx;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.annotation.Import;

@SpringBootApplication
@Import(BeansConfiguration.class)
public class ApiDxApplication {

	public static void main(String[] args) {
		SpringApplication.run(ApiDxApplication.class, args);
	}
}
