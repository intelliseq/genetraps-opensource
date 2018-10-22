package pl.intelliseq.genetraps.api.dx.controllers;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ObjectNode;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.oauth2.common.OAuth2AccessToken;
import org.springframework.security.oauth2.provider.OAuth2Authentication;
import org.springframework.security.oauth2.provider.authentication.OAuth2AuthenticationDetails;
import org.springframework.security.oauth2.provider.token.TokenStore;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.web.bind.annotation.CrossOrigin;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

@RestController
@Log4j2
public class UsersController {

    @Autowired
    private TokenStore tokenStore;

    @Autowired
    private AuroraDBManager auroraDBManager;

    //@CrossOrigin
    @RequestMapping(value = "/user", method = RequestMethod.GET)
    public String user(OAuth2Authentication auth) throws JsonProcessingException {
        //TODO: Choose better option
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        log.debug(String.format("{\"userId\": \"%s\"}", userId));

//        OAuth2AuthenticationDetails details = (OAuth2AuthenticationDetails) auth.getDetails();
//        OAuth2AccessToken accessToken = tokenStore.readAccessToken(details.getTokenValue());
//        log.info(accessToken.getAdditionalInformation().get("user_name").toString());

        String clientId = auth.getOAuth2Request().getClientId();
        log.debug(String.format("{\"client\": \"%s\"}", clientId));

        return new ObjectMapper().writeValueAsString(auroraDBManager.getUserDetails(userId));
    }

    @RequestMapping(value = "groups", method = RequestMethod.GET)
    public String getUsersGroups(OAuth2Authentication auth) {
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        return auroraDBManager.getUsersGroups(userId).toString();
    }
}
