package pl.intelliseq.genetraps.api.security.services;

import com.mongodb.MongoClient;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import com.mongodb.client.model.Filters;
import org.bson.Document;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.core.authority.SimpleGrantedAuthority;
import org.springframework.security.core.userdetails.User;
import org.springframework.security.core.userdetails.UserDetails;
import org.springframework.security.core.userdetails.UserDetailsService;
import org.springframework.security.core.userdetails.UsernameNotFoundException;
import org.springframework.stereotype.Component;

import java.util.Arrays;
import java.util.List;
@Component
public class MongoUserDetailsService implements UserDetailsService{
    @Autowired
    private MongoClient mongoClient;
    @Override
    public UserDetails loadUserByUsername(String username) throws UsernameNotFoundException {

        MongoDatabase database = mongoClient.getDatabase("test");
        MongoCollection<Document> collection = database.getCollection("users");

        Document document = collection.find(Filters.eq("username",username)).first();
        if(document!=null) {
            String password = document.getString("password");
            List<SimpleGrantedAuthority> authorities = Arrays.asList(new SimpleGrantedAuthority("user"));
            return new User(username,password,authorities);
        } else {
            throw new UsernameNotFoundException("User not found");
        }
    }
}
