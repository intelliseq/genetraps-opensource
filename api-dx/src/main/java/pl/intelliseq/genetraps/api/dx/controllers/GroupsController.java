package pl.intelliseq.genetraps.api.dx.controllers;

import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.oauth2.provider.OAuth2Authentication;
import org.springframework.security.oauth2.provider.token.TokenStore;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.RestController;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;
import springfox.documentation.annotations.ApiIgnore;

@RestController
@Log4j2
public class GroupsController {

    @Autowired
    private TokenStore tokenStore;

    @Autowired
    private AuroraDBManager auroraDBManager;

    @RequestMapping(value = "groups", method = RequestMethod.GET)
    public String getGroups(@ApiIgnore OAuth2Authentication auth) {
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        return auroraDBManager.getUsersGroups(userId).toString();
    }

    @RequestMapping(value = "groups", method = RequestMethod.POST)
    public String createGroup(@ApiIgnore OAuth2Authentication auth,
                              @RequestParam String groupName,
                              @RequestParam(required = false, defaultValue = "false") Boolean root) {
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        return null;
    }

    @RequestMapping(value = "group/{id}", method = RequestMethod.GET)
    public String getGroupInfo(@ApiIgnore OAuth2Authentication auth) {
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());

        return null;
    }


}
