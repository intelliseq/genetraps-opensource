package pl.intelliseq.genetraps.api.dx;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.Getter;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public enum UserTest {
    ADMIN(1), DEVIL(3), PSYDUCK(8);

    private JsonNode token;

    @Getter
    private Integer id;

    UserTest(Integer id){
        this.id = id;

        Resource resource = new ClassPathResource(String.format("tokens/%s.json", toString().toLowerCase()));
        try {
            BufferedReader fileReader = new BufferedReader(new FileReader(resource.getFile()));
            ObjectMapper mapper = new ObjectMapper();
            token = mapper.readValue(fileReader, JsonNode.class);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public String getAccessToken(){
        return token.get("access_token").toString();
    }

    public String getUsername(){
        return this.toString().toLowerCase();
    }

}
