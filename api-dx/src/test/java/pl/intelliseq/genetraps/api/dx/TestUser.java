package pl.intelliseq.genetraps.api.dx;

import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

import java.io.FileReader;
import java.io.IOException;

public enum TestUser {
    ADMIN, DEVIL, PSYDUCK;

    private JSONObject token;

    TestUser(){
        Resource resource = new ClassPathResource(String.format("tokens/%s.json", toString().toLowerCase()));

        JSONParser parser = new JSONParser();
        try {
            token = (JSONObject) parser.parse(new FileReader(resource.getFile()));
        } catch (IOException | ParseException e) {
            throw new RuntimeException(e);
        }
    }

    public String getAccessToken(){
        return (String) token.get("access_token");
    }

}
