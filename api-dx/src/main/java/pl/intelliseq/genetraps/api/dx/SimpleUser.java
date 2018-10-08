package pl.intelliseq.genetraps.api.dx;

import lombok.Getter;

public class SimpleUser {
    @Getter
    Integer id;

    @Getter
    String username;

    @Getter
    Boolean root;

    public SimpleUser(Integer id, String username, Boolean root) {
        this.id = id;
        this.username = username;
        this.root = root;
    }

    @Override
    public String toString() {
        return "SimpleUser{" +
                "id=" + id +
                ", username='" + username + '\'' +
                ", root=" + root +
                '}';
    }
}
