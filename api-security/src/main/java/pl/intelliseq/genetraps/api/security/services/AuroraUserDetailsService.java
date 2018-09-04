package pl.intelliseq.genetraps.api.security.services;

import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.core.userdetails.UserDetails;
import org.springframework.security.core.userdetails.UserDetailsService;
import org.springframework.stereotype.Component;
import pl.intelliseq.genetraps.api.security.AuroraDBManager;

@Log4j2
@Component
public class AuroraUserDetailsService implements UserDetailsService{
    @Autowired
    private AuroraDBManager auroraDBManager;
    @Override
    public UserDetails loadUserByUsername(String username) {
        log.info("Getting user details for "+username);
        return auroraDBManager.getUser(username);
    }
}
