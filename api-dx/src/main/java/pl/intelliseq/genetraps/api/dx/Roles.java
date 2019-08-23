package pl.intelliseq.genetraps.api.dx;

import lombok.Getter;

public enum Roles {
    READER(1), WRITER(2), ADMIN(3);

    @Getter
    Integer id;

    Roles(Integer id) {
        this.id = id;
    }
}
