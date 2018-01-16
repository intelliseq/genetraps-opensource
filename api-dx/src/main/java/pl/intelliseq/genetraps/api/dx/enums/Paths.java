package pl.intelliseq.genetraps.api.dx.enums;

public enum Paths {

    RUN_URL_FETCHER(System.getProperty("user.dir") + "/src/main/resources/bash/runUrlFetcher.bash"),
    GET_URL_FETCHER_STATE(System.getProperty("user.dir") + "/src/main/resources/bash/getUrlFetcherState.bash");

    private String path;

    private Paths(String path) {this.path = path; }

    public String getPath() {return this.path; }

}
