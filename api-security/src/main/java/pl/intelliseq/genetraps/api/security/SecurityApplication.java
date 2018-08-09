package pl.intelliseq.genetraps.api.security;

import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.builder.SpringApplicationBuilder;

@SpringBootApplication
public class SecurityApplication {

	public static void main(String[] args) {
		new SpringApplicationBuilder(SecurityApplication.class)
				.properties("spring.config.name=authserver")
                .run(args);
	}
}
