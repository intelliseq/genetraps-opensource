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
import org.springframework.web.bind.annotation.*;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;

@RestController
@Log4j2
public class GroupsController {

    @Autowired
    private TokenStore tokenStore;

    @Autowired
    private AuroraDBManager auroraDBManager;

    @RequestMapping(value = "groups", method = RequestMethod.GET)
    public String getGroups(OAuth2Authentication auth) {
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        return auroraDBManager.getUsersGroups(userId).toString();
    }

    @RequestMapping(value = "groups", method = RequestMethod.POST)
    public String createGroup(OAuth2Authentication auth,
                              @RequestParam String groupName,
                              @RequestParam(required = false, defaultValue = "false") Boolean root){
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        return null;
    }

    @RequestMapping(value = "group/{id}", method = RequestMethod.GET)
    public String getGroupInfo(OAuth2Authentication auth){
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        return null;
    }


}
