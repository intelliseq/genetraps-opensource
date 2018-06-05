package pl.intelliseq.explorare.model.hpo;

import java.util.List;

public class ObsoleteHpoTerm {

    private String id;
    private String name;
    private List<String> synonyms;

    public ObsoleteHpoTerm(String id, String name, List<String> synonyms) {
        this.id = id;
        this.name = name;
        this.synonyms = synonyms;
    }

    public String getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public List<String> getSynonyms() {
        return synonyms;
    }

}
