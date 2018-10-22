package pl.intelliseq.genetraps.api.dx.controllers;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ObjectNode;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.oauth2.provider.OAuth2Authentication;
import org.springframework.web.bind.annotation.*;
import pl.intelliseq.genetraps.api.dx.Roles;
import pl.intelliseq.genetraps.api.dx.SimpleUser;
import pl.intelliseq.genetraps.api.dx.exceptions.ForbiddenException;
import pl.intelliseq.genetraps.api.dx.helpers.AuroraDBManager;
import pl.intelliseq.genetraps.api.dx.helpers.DxApiProcessManager;
import pl.intelliseq.genetraps.api.dx.helpers.FilesManager;

@RestController
@Log4j2
public class PrivilegesController {
    @Autowired
    private DxApiProcessManager processManager;

    @Autowired
    private FilesManager filesManager;

    @Autowired
    private AuroraDBManager auroraDBManager;

    @RequestMapping(value = "/user/privileges", method = RequestMethod.GET)
    public String getUserPrivileges(OAuth2Authentication auth){
        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());
        ObjectNode result = new ObjectMapper().createObjectNode();

        //map to json
        auroraDBManager.getUserPrivileges(userId).forEach((key, value) -> result.put(key.toString(), value.toString()));
        return result.toString();
    }

    @RequestMapping(value = "user/privileges", method = RequestMethod.POST)
    public String giveUserPriviliges(
            OAuth2Authentication auth,
            @RequestParam Integer targetUserId,
            @RequestParam Integer sampleId,
            @RequestParam Roles role){

        Integer userId = Integer.valueOf(auth.getUserAuthentication().getPrincipal().toString());
        SimpleUser user = auroraDBManager.createNewSimpleUser(userId);
        Roles loggedUserRole = auroraDBManager.getUserPrivilegesToSample(user.getId(), sampleId);

        if (!loggedUserRole.equals(Roles.ADMIN)){
            throw new ForbiddenException("User don't have permissions to share sample");
        }

        return String.format("{\"response\":%s}", auroraDBManager.setUserPrivilegesToSample(targetUserId, sampleId, role));
    }

}
