package pl.intelliseq.genetraps.api.dx;

import lombok.Getter;
import lombok.Setter;

@Getter
public class User {
    @Setter
    private Integer id;
    private final String lastName;
    private final String firstName;
    private final String email;

    private final String userName;
    private final Boolean root;

    public static Builder builder(){
        return new Builder();
    }

    private User(Builder builder){
        this.id = builder.id;
        this.lastName = builder.lastName;
        this.firstName = builder.firstName;
        this.email = builder.email;
        this.userName = builder.userName;
        this.root = builder.root;
    }

    @Override
    public String toString() {
        return "User{" +
                "id=" + id +
                ", lastName='" + lastName + '\'' +
                ", firstName='" + firstName + '\'' +
                ", email='" + email + '\'' +
                ", userName='" + userName + '\'' +
                ", root=" + root +
                '}';
    }

    public static class Builder{
        private Integer id;
        private String lastName;
        private String firstName;
        private String email;

        private String userName;
        private Boolean root = false;

        public User build(){
            return new User(this);
        }

        public Builder withId(Integer id){
            this.id = id;
            return this;
        }

        public Builder withLastName(String lastName){
            this.lastName = lastName;
            return this;
        }

        public Builder withFirstName(String firstName){
            this.firstName = firstName;
            return this;
        }

        public Builder withEmail(String email){
            this.email = email;
            return this;
        }

        public Builder withUserName(String userName){
            this.userName = userName;
            return this;
        }

        public Builder withRoot(Boolean root){
            this.root = root;
            return this;
        }
    }
}
